#ifndef EE_H
#define EE_H
#include"CGraph.h"
#include <ilcplex/ilocplex.h>


// 4+x^2
double EEdictor(CGraph *G,vector<demand> & req,int ornum,double OPEN){
	IloEnv env;
	IloModel model(env);
	IloCplex EEsolver(model);

	int num = req.size();	
	IloArray<IloIntVarArray> x(env, num); 
	for(int d = 0; d < num; d++)
		x[d] = IloIntVarArray(env, G->m, 0, 1); 	

	// 对每个点进行流量守恒约束  
	for(int d = 0; d < num; d++){
		for(int i = 0; i < G->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < G->adjL[i].size(); k++) // 出度边
				constraint += x[d][G->adjL[i][k]->id];
			for(unsigned int k = 0; k < G->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][G->adjRL[i][k]->id];

			if(i == req[d].org)
				model.add(constraint == 1);
			else if(i == req[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}
	}

	for(int i = 0; i < G->m; i++){
		IloExpr constraint(env);
		for(int d = 0; d <  num; d++)
			constraint += req[d].flow*x[d][i];
		model.add( constraint <= G->Link[i]->capacity );  
	}

	//优化目标 节能 OPEN+X^2
	IloExpr cost(env);
	for(int i = 0; i < G->m; i++){
		IloExpr load(env);
		IloIntVarArray tmp(env,num,0,1);
		for(int d = 0; d < num; d++){
			load += req[d].flow*x[d][i];
			tmp[d] = x[d][i];
		}
		//cost += ( IloPower(load,2) + IloMax(tmp)*OPEN );
		cost += ( load*load + IloMax(tmp)*OPEN );
	}
	model.add(IloMinimize(env,cost));

	EEsolver.setOut(env.getNullStream());
	double obj = INF;
	if(EEsolver.solve()){
		obj = EEsolver.getObjValue(); //energy
		//delay
		for(int i=0;i<G->m;i++){  
			double loadc = 0;
			for(int d=0;d < num;d++)
				loadc += EEsolver.getValue(x[d][i])*req[d].flow;
			G->Link[i]->latency = linearCal(loadc,G->Link[i]->capacity);
		}
		double del = 0;
		for(int d=0;d < ornum;d++){
			for(int i=0;i<G->m;i++)
				del +=  EEsolver.getValue(x[d][i])*G->Link[i]->latency;
		}
		G->delay = del;
		cout << "EE\t"<< obj <<"\t"<<del<<endl;
	}
	for(int i = 0; i < req.size(); i++)
		x[i].end();
	x.end();
	env.end();
	return obj;
}

double ORdictor(CGraph *G,vector<demand> eq,int ornum,double OPEN){    
	IloEnv env;
	IloModel model(env);
	IloCplex ORsolver(model);

	int num = eq.size();
	IloArray<IloNumVarArray> x(env, num); 
	for(int d = 0; d < num; d++)
		x[d] = IloNumVarArray(env, G->m, 0, 1); // num * G->m 的二维矩阵


	if(!CONSTANT){
		IloNumVarArray cost(env,G->m,0.0 ,IloInfinity); 
		IloExpr res(env);
		for(int i=0;i<G->m;i++)
		{
			IloNumExpr load(env);
			for(int d=0;d < ornum;d++)
				load+=x[d][i]*eq[d].flow;
			double c=G->Link[i]->capacity;
			model.add(cost[i] >= load);
			model.add(cost[i] >= 3*load-2*c/3);
			model.add(cost[i] >= 10*load-16*c/3);
			model.add(cost[i] >= 70*load-178*c/3);
			res += cost[i];
		}
		model.add(IloMinimize(env, res));

		//for(int i=0;i<G->m;i++){
		//	IloExpr load(env);
		//	for(int d=0;d<ornum;d++){ 
		//		load += eq[d].flow*x[d][i];
		//	}
		//	model.add(cost[i]>=load*load/(G->Link[i]->capacity));   // [0,1/3]
		//	model.add(cost[i]>=3.0*load*load/(G->Link[i]->capacity)-(2.0/3.0));  //[1/3,2/3]
		//	model.add(cost[i]>=10.0*load*load/(G->Link[i]->capacity)-(16.0/3.0)); // [2/3,9/10]
		//	model.add(cost[i]>=70.0*load*load/(G->Link[i]->capacity)-(178.0/3.0));  //[9/10,1]
		//	model.add(cost[i]>=500*load*load/G->Link[i]->capacity-(double)(1468/3)); // [1,11/10]
		//	model.add(cost[i]>=5000*load*load/G->Link[i]->capacity-(double)(19468/3)); //[11/10,...]
		//   res += cost[i];	
		//}
		//model.add(IloMinimize(env,res));
	}
	else{
		IloExpr delay(env);
		for(int ij = 0; ij  < G->m; ij++){		
			IloNumExpr ep1(env);
			for(int d = 0; d < ornum; d++){
				ep1 += x[d][ij]*G->Link[ij]->dist;  //d之和			
			}
			delay += ep1;		 
		}
		model.add(IloMinimize(env, delay)); 
	}

	//约束  流量守恒约束  对每个点进行流量守恒约束  
	for(int d = 0; d < num; d++)
		for(int i = 0; i < G->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < G->adjL[i].size(); k++) // 出度边
				constraint += x[d][G->adjL[i][k]->id];
			for(unsigned int k = 0; k < G->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][G->adjRL[i][k]->id];
			// 出 - 入
			if(i == eq[d].org)
				model.add(constraint == 1);
			else if(i == eq[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}

		//带宽约束
		for(int i = 0; i < G->m; i++){
			IloExpr constraint(env);
			for(int d = 0; d <  num; d++)
				constraint += eq[d].flow*x[d][i];
			model.add(constraint <= G->Link[i]->capacity);  
		}

		ORsolver.setOut(env.getNullStream());
		double obj = INF;
		ORsolver.solve();
		if(ORsolver.getStatus() == IloAlgorithm::Infeasible)
			env.out() << "OR Dictator No Solution" << endl;
		else{
			obj = ORsolver.getObjValue();
			// 计算 energy
			double util = 0;
			for(int i = 0; i < G->m; i++){  
				double load = 0;
				double one = 0;
				for(int d = 0; d < num; d++){
					load += ORsolver.getValue(x[d][i])*eq[d].flow;
					one = max(one,ORsolver.getValue(x[d][i]));
				}
				util += (one*OPEN + load*load);
				//util += (one*OPEN + load);
			}
			G->energy = util;

			// delay
			double del = 0;
			for(int i=0;i<G->m;i++){  
				double loadc = 0;
				for(int d = 0;d < ornum;d++)
					loadc += ORsolver.getValue(x[d][i])*eq[d].flow;
				del += linearCal(loadc,G->Link[i]->capacity);
			}
			obj = del;
			cout << "OR\t"<< util<<"\t"<<obj <<endl;
		}
		for(int i = 0; i < eq.size(); i++)
			x[i].end();
		x.end();
		env.end();
		return obj;
}

double CGraph::EE(int id,int s,int t,double dm,bool needpath,double OPEN){
	vector<int> p, flag;
	vector<double> d;//energy
	for(int i = 0; i < n; i++){
		p.push_back(-1);
		flag.push_back(0);
		d.push_back(INF);
	}

	d[s] = 0;
	int cur = s;
	do{	
		flag[cur] = 1;
		for(unsigned int i = 0; i < adjL[cur].size(); i++){
			CEdge *e = adjL[cur][i];
			//e->latency = pow((e->use+dm),2)-pow(e->use,2);
			e->latency = pow((e->use+dm),2);
			if( 0 == e->use )
				e->latency += OPEN;
			if(e->capacity - e->use >= dm && d[e->head] > d[e->tail] + e->latency){
				d[e->head] = d[e->tail] + e->latency;
				p[e->head] = e->id;
			}
		}
		cur = -1;
		for(int i = 0; i < n; i++)
			if(!flag[i] && (cur == -1 || d[cur] > d[i] ))
				cur = i;
	}while(cur != -1);

	cur = t;
	do{
		if(p[cur] == -1)
			break;
		Link[p[cur]]->use += dm;
		cur = Link[p[cur]]->tail;
	}while(cur != s);

	if(needpath){
		reqPathID[id].clear();
		cur = t;
		do{
			if(p[cur] == -1)
				break;
			this->reqPathID[id].push_back(p[cur]);
			cur = Link[p[cur]]->tail;
		}while(cur != s);
		reverse(reqPathID[id].begin(),reqPathID[id].end());
	}

	return d[t];
}

#endif