#ifndef EVOLUTIONBIT_H
#define EVOLUTIONBIT_H
#include"DFS.h"
#include"EE.h"

// 个体 individual
class evoluDivbit{
private:
	static const int MUT = 6; //变异位数
	static const int HER = 15; //学习位数
	double consider;
	CGraph *G;
	CGraph *GOR;
	vector<demand> *dem; 
	vector<demand> *demor; 
public:
	vector<vector<int>> x;//解  各个req经过的路径编号 二进制编码 req*4
	double ability; //个体能力
	double delay, energy;//个体能力值的体现
	double delaybase;
	double energybase;
	double OPEN;
	// 构造函数
	evoluDivbit() {;}
	double GAabilityNoCplex();
	void Caldelay();
	void Init(){
		ability = MIN;
		delay = INF;
		energy = INF;
	}
	//m为req的数量
	evoluDivbit(int m, CGraph *g, CGraph *gor,vector<demand> *d,vector<demand> *dor,double energybest,double delaybest, double con,double open){
		Init();
		x.resize(m);//
		G = g;
		GOR=gor;
		dem = d;
		demor = dor;
		delaybase=delaybest;
		energybase = energybest;
		OPEN = open;
		consider = con;
		randNature(); 
	}

	evoluDivbit(vector<vector<int> > &tx, CGraph *g,CGraph *gor, vector<demand> *d, vector<demand> *dor,double energybest,double delaybest,double con,double open){
		Init();
		x.clear();
		G = g;
		GOR=gor;
		dem = d;
		demor=dor;
		delaybase=delaybest;
		energybase = energybest;
			OPEN = open;
		consider = con;
		for(unsigned int i = 0; i < tx.size(); i++)
			x.push_back(tx[i]);
		// 也可以用  x = tx ; 代替上面的操作
	}

	/////////////交配杂交  crossover
	evoluDivbit crossover(evoluDivbit other){
		vector<vector<int> > nx;	
		/*
		//分成两截互换 在num处截断  即是单点杂交
		int n=dem->size();
		int num=rand()%(n/2);
		for(int i=0;i<num;i++) 
		nx.push_back(x[i]);
		for(unsigned int ij=num;ij<x.size();ij++)
		nx.push_back(other.x[ij]);	
		return evoluDivbit(nx, G,GOR, dem,demor);	
		*/			
		vector<vector<int> > onezero;   //多点交叉  随机产生一串0-1序列  为0不交叉  为1交叉
		for(unsigned int i=0;i<x.size();i++){
			vector<int> bit;
			for(int j=0;j<4;j++)  
				bit.push_back(rand()%2);
			onezero.push_back(bit); 
		}
		for(unsigned int n=0;n<x.size();n++)
			nx.push_back(x[n]);
		for(unsigned int i=0;i<onezero.size();i++)
			for(int j=0;j<4;j++){  
				if(onezero[i][j]==1)
					nx[i][j]=other.x[i][j];
			}
			return evoluDivbit(nx, G,GOR, dem,demor,other.energybase,other.delaybase,other.consider,other.OPEN);	
	}

	//将二进制转为十进制  路径编号的具体值
	int Decoding(vector<int> &bit){
		int z=0;
		z = bit[0]*1+bit[1]*2+bit[2]*4+bit[2]*8;//// 16条路
		//z = bit[0]*1+bit[1]*2+bit[2]*4;//// 8条路
		return z;
	}

	// *********  计算能力  *********//
	void calAbility(){
		G->reqPathNo.clear(); //清空上一个个体解的路径编号
		for(unsigned int i=0;i<x.size();i++) 
			G->reqPathNo.push_back(Decoding(x[i])); 
		ability = GAabilityNoCplex();	
	}

	void randNature(){
		for(unsigned int i = 0; i < x.size(); i++){
			vector<int> bit;
			for(int b=0;b<4;b++)        ///////////////////////   16 ： 4      8 ： 
				bit.push_back(rand()%2);// 4位编码 随机产生0或1  表示路径编号
			x[i] = bit; 
		}
	}

	void mutation(){	//变异
		int i=0;
		while(i<MUT){
			int row= rand()% dem->size();//产生需要变异的req路径编号
			int col=rand()%4;//产生需要变异的位置   ///////////////////////   16 ： 4      8 ： 3
			if(x[row][col]==0) x[row][col]=1;
			else if(x[row][col]==1) x[row][col]=0; 
			i++;
		}
	}

	/////////学习历史最优解
	void culture(evoluDivbit hero){
		int n=0;
		while(n<HER){
			int j=rand()%x.size();
			for(unsigned int i = 0; i < x.size(); i++){
				if(i==j) x[i] =  hero.x[i];
			}
			n++;
		}	
	}

};

bool Cmp2(evoluDivbit a, evoluDivbit b){
	return a.ability > b.ability; //ability和适应度成反比，ability越小（最小化问题），适应度越大
}

//种群 种群由个体组成
class evoluPopubit{
private:
	static const int YEAR = 70;
	static const int YEARCUL = 50;
	static const int NOHERO = 30;
	double pm;
	vector<evoluDivbit> popu;
	CGraph *G;
	CGraph *GOR;
	vector<demand> *dem;
	vector<demand> *demor;
	FILE *herofile;
public:
	evoluDivbit hero;
	// n 种群大小  m:req数目 初始化种群
	evoluPopubit(int n, int m, CGraph *g, CGraph *gor,vector<demand> *d,vector<demand> *dor,double energybest,double delaybest,double con,double OPEN){
		popu.clear();
		pm = 0.25;
		G = g;
		GOR = gor;
		dem = d;
		demor = dor;
		for(int i = 0; i < n; i++){
			evoluDivbit divi(m, G, GOR,dem,demor,energybest,delaybest,con,OPEN);
			popu.push_back(divi);
		}
		hero = evoluDivbit(m, G,GOR, dem,demor,energybest,delaybest,con,OPEN);
		herofile = fopen("outputFile//hero.txt","a");
	}
	///////////  轮盘赌和
	int wheelchose(double sum){
		double m=0;
		double r=rand()%1000*0.001;
		int res = rand()%popu.size();
		for(unsigned int i=0;i<popu.size();i++){
			m += popu[i].ability/sum;
			if(r<=m) {
				//	cout << i << "   ";
				return i;
			}
		}
		return res; //如果for循环中没有返回   // 解决 ：“evoluPopubit::wheelchose”: 不是所有的控件路径都返回
	}

	/////////// 种群进化 //////////////////
	evoluDivbit evolution(){
		int nohero=0;//记录进化情况
		fprintf(herofile,"Start:\n ");
		vector<int> h0;
		for(int b=0;b<4;b++)  
			h0.push_back(0);
		for(unsigned int i = 0; i < hero.x.size(); i++)
			hero.x[i] = h0;//开始的hero设为每个req都走各自路径集合的第0条路
		hero.calAbility();
		if(hero.ability>MIN) fprintf(herofile, "%f\t%f\t%f\n", hero.energy,hero.delay,hero.ability);
		//评价每个个体的能力
		for(unsigned int i = 0; i < popu.size(); i++)
			popu[i].calAbility();
		//sort(popu.begin(), popu.end(), Cmp2);//按个体优劣排列以进行强强交配
		//繁殖的代数
		for(int curYear = 1; curYear <= YEAR; curYear++){
			//轮盘赌和选择n个个体
			int n = popu.size(), getMore = 0;
			double sum=0;
			for( int i=0;i<n;i++)
				sum += popu[i].ability;
			vector<evoluDivbit> chosepopu;
			chosepopu.clear();
			//cout<< " ch   " <<endl;
			for(int i=0;i<n;i++) //n种群大小
			{
				int ch=wheelchose(sum);
				// cout<<"  "<<ch<<"  ";
				chosepopu.push_back(popu[ch]);
			}
			popu.clear();
			for(unsigned int i=0;i<chosepopu.size();i++)
				popu.push_back(chosepopu[i]);
			vector<evoluDivbit> sons;//子孙
			sons.clear();
			//杂交
			for(unsigned int i = 0; i+1 < popu.size(); i+=2){
				sons.push_back(popu[i].crossover(popu[i+1]));
				sons.push_back(popu[i+1].crossover(popu[i]));
				sons.push_back(popu[i].crossover(popu[i+1]));
			}
			// sons.push_back(hero);  //精英保留策略
			int m = sons.size();//子孙的个数
			for(int i = 0; i < m; i++){
				double p=rand()%100*0.01;
				if(p<pm) 
					sons[i].mutation();//个体变异以概率pm
				if(curYear > YEARCUL)
					sons[i].culture(hero); //向hero靠近
				sons[i].calAbility();
			}
			sort(sons.begin(), sons.end(), Cmp2);
			popu.clear(); //清空上一代种群
			//printf("Year %d\n ", curYear);
			for(int i = 0; i < n; i++){
				//printf("%.6f ", sons[i].ability);
				popu.push_back(sons[i]); //保持种群数量不变为n,将能力排在前n的放入下一代种群
				if(abs(sons[i].ability - hero.ability) < 1e-4 ){	//判断相等
					//cout<< sons[i].ability <<" ====== "<< hero.ability;
					continue;
				}
				else if(sons[i].ability > hero.ability){
					//cout<< sons[i].ability <<"  >>>>>>  "<< hero.ability;
					hero = sons[i]; //选出能力最好的作为hero
					getMore = 1;
				}
			}
			if(getMore){
				fprintf(herofile, "Year %d: find hero \n", curYear);
				fprintf(herofile,"%f\t%f\t%f\n", hero.energy,hero.delay,hero.ability);
			}
			else nohero++;
			if(nohero> NOHERO){
				break;
			}
		}
		fprintf(herofile,"end\n\n\n\n");
		fclose(herofile);
		return hero;
	}
};

double evoluDivbit::GAabilityNoCplex(){ 
	this->G->clearOcc(); 
	for(unsigned int d=0;d<G->reqPathNo.size();d++){
		int num=G->reqPathNo[d]; 
		//vector<CEdge*> lsEg=G->reqlistPath[d][num]->listEdge; //KSP
		vector<CEdge*> lsEg=G->reqlistPath[d][num]; //DFS
		vector<CEdge*>::iterator it,end=lsEg.end();
		for(it=lsEg.begin();it!=end;it++){
			(*it)->use += (*dem)[d].flow;
		}
	}

	this->GOR->clearOcc();
	Caldelay();
	/*for (int i = 0; i < GOR->m; i++)
		cout << GOR->Link[i]->latency<<" ";
	cout <<endl<<endl;*/
	double del = 0;
	for (unsigned int d = 0; d < demor->size(); d++){		
			del += GOR->dijkstraOR((*demor)[d].org, (*demor)[d].des, (*demor)[d].flow);
	}
	//cout <<endl;

	//新的流量矩阵计算TE的目标
	vector<demand> req;
	for (int i = 0; i < GOR->m; i++){
		if( GOR->Link[i]->use > 0 )
			req.push_back(demand(GOR->Link[i]->tail, GOR->Link[i]->head, GOR->Link[i]->use));
	}
	int ornum = req.size();
	for (unsigned int i = demor->size(); i < dem->size(); i++){
		req.push_back((*dem)[i]);
	}

	// 计算EE的目标值
	G->clearOcc(); 
	double ee,del2;
	heuristicEE(G,req,ornum,ee,del2,OPEN);
	
	//cout << del <<" "<<del2<<" ***;   ";
	double ab = MIN;
	if ( del+1e-5 <= INF && ee+1e-5 <= INF){
		this->energy = ee;
		this->delay = del;	
		ab = this->energybase/ee + this->delaybase*this->consider / del;  
	}
	return ab;
}


void evoluDivbit::Caldelay(){
	for(int ij=0;ij<G->m;ij++){
		G->Link[ij]->latency = linearCal(G->Link[ij]->use,G->Link[ij]->capacity);
		//cout << G->Link[ij]->latency<<"  ";
	}
	//cout <<endl;
	int totalnum = (*dem).size();
	for(int m=0;m<GOR->m;m++){	
		bool flag=false;
		double dmin=0;
		for(int d=0;d<totalnum;d++){
			if (GOR->Link[m]->tail == (*dem)[d].org && GOR->Link[m]->head == (*dem)[d].des && ((*dem)[d].flow>0)){
				if(flag==false){
					flag=true;
					//vector<CEdge*> lg=G->reqlistPath[d][G->reqPathNo[d]]->listEdge; //KSP
					vector<CEdge*> lg=G->reqlistPath[d][G->reqPathNo[d]]; //DFS
					for(vector<CEdge*>::iterator lgi=lg.begin();lgi!=lg.end();lgi++) 
							dmin += (*lgi)->latency;
				}
			}
			if(flag==true){
				GOR->Link[m]->latency = dmin;	
				break;
			}
		} 
		if(flag == false)
			GOR->Link[m]->latency = G->dijkstraOR(GOR->Link[m]->tail,GOR->Link[m]->head,0);
	} 
}

#endif
