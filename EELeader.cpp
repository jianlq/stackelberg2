#include"evolutionbit.h"
#include"EE.h"
#include"nash.h"

// 在友好系数0.4下，对于不同的启动能耗，测试nash，stackelberg，EE，OR
// 每个启动能耗下进行30个测试案例的测试，每个启动能耗最终的结果取30次的平均
int main(){
	srand((unsigned)time(NULL));
	int Time = 3;
	int CASEnum= 30;	
	vector<double> STARTUP;
	STARTUP.push_back(0);
	STARTUP.push_back(10);
	STARTUP.push_back(100);
	STARTUP.push_back(500);
	STARTUP.push_back(1000);
	STARTUP.push_back(2000);
	STARTUP.push_back(4000);
	STARTUP.push_back(10000);

	vector<double>CONSIDER;
	int CON_VALUE = 0;
	double c = 0.4;
	for(int i=0;i <= CON_VALUE;i++){
		CONSIDER.push_back(c);
		c += 0.2;
	}

	int LOOP  = 100;

	for(int i =0;i<Time;i++){
		for(int start = 0; start < STARTUP.size();start++){
			FILE *out = fopen("outputFile//eeor.csv", "a");
			FILE *res = fopen("outputFile//result.csv", "a");
			FILE *nash = fopen("outputFile//nash.csv", "a");

			int conN = CONSIDER.size();
			vector<double> see(conN,0) ;
			vector<double> sor(conN,0);

			vector<double> dicee(conN,0) ;
			vector<double> diceeor(conN,0) ;

			vector<double> dicoree(conN,0);
			vector<double> dicor(conN,0);

			vector<double> nashee(conN,0) ;
			vector<double> nashor(conN,0);

			vector<int> successCase (conN, 0) ;
			vector<int> flag(conN,1);

			for(int casenum = 0; casenum < CASEnum; casenum++){
				genGraph(9,50,"inputFile//graph.txt");
				genGraphOR(9,4,9,"inputFile//graphOR.txt");
				CGraph *G = new CGraph("inputFile//graph.txt");
				CGraph *GOR = new CGraph("inputFile//graphOR.txt");

				vector<demand> eqOR;
				vector<demand> eqTE;
				vector<demand> eqbase;

				//eqbase.clear();//background流 
				for(int i = 0; i < BGNUM; i++){
					int s = rand()%G->n, t;
					do{
						t = rand()%G->n;
					}while( s == t || G->canNotReach(s,t));
					eqbase.push_back(demand(s, t, rand()%(MAXDEMAND)+1));
				}

				////Overlay  产生demand
				//eqOR.clear(); 
				for(int i =0 ;i<GOR->m;i++)
					if(G->canNotReach(GOR->Link[i]->tail, GOR->Link[i]->head))
						continue;
					else
						eqOR.push_back(demand(GOR->Link[i]->tail,GOR->Link[i]->head,rand()%(MAXDEMAND)+1));

				//eqTE.clear(); 
				for(unsigned int i=0;i<eqOR.size();i++)
					eqTE.push_back(eqOR[i]);
				int ornum = eqTE.size();
				for(unsigned int i=0;i<eqbase.size();i++)
					eqTE.push_back(eqbase[i]);

				double eedic = 0, ordic = 0;
				G->clearOcc();
				eedic = EEdictor(G,eqTE,ornum,STARTUP[start]);

				G->clearOcc();
				ordic = ORdictor(G,eqTE,ornum,STARTUP[start]);

				if( (eedic + 1e-5 >= INF)  || ( ordic + 1e-5 >=INF) )
					break;

				G->clearOcc();
				if(!G->GAinit(eqTE)){
					cout << "*****GA init failed***** "<<endl;
					break;
				}

				
				for(int con = 0;con < CONSIDER.size();con++){

					int n = 150;//种群个体数目
					int m = eqTE.size();
					evoluPopubit popubit(n,m,G,GOR,&eqTE,&eqOR,eedic,ordic,CONSIDER[con],STARTUP[start]);
					(popubit).evolution();
					cout<<"S\t"<<popubit.hero.energy<<"\t"<<popubit.hero.delay<<endl;

					if( (popubit.hero.energy +1e-5) >= INF ||  (popubit.hero.delay + 1e-5) > INF ){
						flag[con] = 0;
						break;
					}

					//// nash	
					int nacase = 0;
					double loopnashee=0,loopnashor=0;
					fprintf(out,"\n nash \n");
					for(int i =0;i<LOOP;i++){
						G->clearOcc();
						double ee = NashEE(G,GOR,eqTE,STARTUP[start]);
						if(ee + 1e-5 >= INF){
							fprintf(nash,"NashEE unfeasible\n");
							break;
						}
						double del = NashOR(GOR,eqOR);
						if(del + 1e-5 >= INF){
							fprintf(nash,"NashOR unfeasible\n");
							break;
						}
						eqTE.clear();
						for(int i=0;i<GOR->m;i++){
							if(GOR->Link[i]->use>0)
								eqTE.push_back(demand(GOR->Link[i]->tail,GOR->Link[i]->head,GOR->Link[i]->use));
						}
						for(unsigned int i=0;i<eqbase.size();i++)
							eqTE.push_back(eqbase[i]);
						GOR->clearOcc();
						nacase++;
						loopnashee += ee;
						loopnashor += del;
						fprintf(nash,"%f,%f\n",ee,del);
					}
					fclose(nash);
				

					if(flag[con]){	
						fprintf(out,"EE,%f,%f,%f\n",STARTUP[start],eedic,G->delay);
						fprintf(out,"OR,%f,%f,%f\n",STARTUP[start],G->energy,ordic);
						fprintf(out,"S,%f,%f,%f\n",STARTUP[start],popubit.hero.energy,popubit.hero.delay);
						fprintf(out,"Nash,%f,%f,%f\n",STARTUP[start],loopnashee/nacase,loopnashor/nacase);

						successCase[con] += 1;
						dicee[con] += eedic;
						diceeor[con] += G->delay;

						dicoree[con] += G->energy;
						dicor[con] += ordic;						

						see[con] += popubit.hero.energy;
						sor[con] += popubit.hero.delay;

						nashee[con] += (loopnashee/nacase);
						nashor[con] += (loopnashor/nacase);
						fclose(out);
					}
				}
								
				delete G;
				delete GOR;
			}
			fprintf(res, "%f\n",STARTUP[start]);
			for(int con = 0;con < CONSIDER.size();con++){
				fprintf(res, "EE,%f,%f,%f\n",CONSIDER[con],dicee[con]/successCase[con],diceeor[con]/successCase[con]); 
				fprintf(res, "OR,%f,%f,%f\n",CONSIDER[con],dicor[con]/successCase[con],dicoree[con]/successCase[con]); 
				fprintf(res, "S,%f,%f,%f\n",CONSIDER[con],see[con]/successCase[con],sor[con]/successCase[con]); 
				fprintf(res, "nash,%f,%f,%f\n",CONSIDER[con],nashee[con]/successCase[con],nashor[con]/successCase[con]); 
			}
			fclose(res);
		}
	}
	system("pause");
	return 0;	
}