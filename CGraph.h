#ifndef CGRAPH_H
#define CGRAPH_H
#include"Common.h"

class CEdge{
public:
	int  id, tail,head;
	double use;
	double dist,latency,load,capacity;
	CEdge(){ }
	CEdge(int i, int a, int b, double c,double d){
		id = i;
		tail = a;
		head = b;
		dist = c;
		capacity = d;
		load = 0; 
		latency = INF;
		use  = 0;
	}

	~CEdge(){}
	double getWeight(){
		return dist;
	}
	int getHead(){
		return head;
	}
	int getTail(){
		return tail;
	}
	double getCap(){
		return capacity;
	}
	bool operator<(CEdge& x)//运算符重载
	{
		if(dist < x.dist)
			return 1;
		else
			return 0;
	}
};

class CVertex{
    public:
	double d;
	int p;
	int ID;
	CVertex(){d = INF; p = NULL;}
	CVertex(int i){ID=i;d = INF; p = NULL;}
	~CVertex();
};
bool pVertexComp ( CVertex* x, CVertex* y )
{    if ( x->d < y->d ) return 1;
	return 0;
};
class Status{
    public:
    int ID;
    double d;
    list<int> passby;
    Status* parent;

    Status(){ID = NULL; d = INF ; };
    Status(int id, double w){ID = id; d = w ; passby.clear();};
    Status(int id, double w, Status* p){ID = id; d = w ; parent = p; passby.clear();};
    Status(int id, double w, Status* p, list<int> listVID){ID = id; d = w ; parent = p; passby = listVID;};

    ~Status(){;};
};
bool pStatusComp ( Status* x, Status* y )
{    if ( x->d < y->d ) return 1;
	return 0;
};


class CPath{
    public:
    double length;
    vector<CEdge*> listEdge;
    list<CVertex*> listVertex;

    CPath(){length = 0; listEdge.clear(); listVertex.clear();};
    CPath(Status *beginning, Status *ending, map<int, CVertex*> mapVID_Vertex, map<int, vector<CEdge*> > mapVID_listEdge);
    ~CPath(){;};
};
bool pPathComp ( CPath* x, CPath* y )
{    if ( x->length < y->length ) return 1;
	return 0;
};
CPath::CPath(Status *beginning, Status *ending, map<int, CVertex*> mapVID_Vertex,map<int, vector<CEdge*> > mapVID_listEdge)
{
    length = ending->d;
    Status *status = ending;
    while(1)
    {
        int id = status->ID;
        //cout<<id<<"  ";
        listVertex.push_front(mapVID_Vertex[id]);
        if(status == beginning)
            break;
        status = status->parent;
       vector<CEdge*>::iterator it;
        for(it = mapVID_listEdge[status->ID].begin(); it != mapVID_listEdge[status->ID].end(); it++)
        {
            if((**it).getHead() == id)
                listEdge.push_back(*it);
        }
    }
}


class demand{
public:
	int org;
	int des;
	double flow;
	demand(int a,int b ,double c)
	{
		org=a;
		des=b;
		flow=c;
	}
};

class CGraph{
private:
	void dfs(int cur){
		visit[cur] = 1;
		for(unsigned int i = 0; i < adjL[cur].size(); i++)
			if(!visit[adjL[cur][i]->head])
				dfs(adjL[cur][i]->head);
	}
public:
	int n, m;
	vector<CEdge*> Link;
	vector<int> ver;	// 所有的顶点
	vector<vector<CEdge*> > adjL, adjRL; //出、入度

	// ksp
	map<int, CVertex*> mapVID_Vertex;	
	map<int, vector<CEdge*> > mapVID_listEdge; // 记录与顶点关联的出度边	
	list<Status*> listTemp;///记录临时状态指针	
	list<Status*> listSolu;///记录被标记为解的状态	
	vector<CPath*> listPath; ///记录path集合
	void KSP(int s, int t, unsigned int k);


	////遗传 
	vector<int> reqPathNo; //所有req走的路径编号
	bool GAinit(vector<demand>& req);//初始化每个req的路径
	void myDFS(int s,int t);
	int DFSflag;
	vector<vector<vector<CEdge*>>> reqlistPath; //DFS
	//vector<vector<CPath*> > reqlistPath;//KSP
	vector<int> pathver;

	// dijkstra
	vector<vector<int>> reqPathID;
	double dijkstra(int id,int s, int t, double dm ,bool ORselfish,bool needpath);
	double dijkstraOR(int s, int t,double dm);
	double EE(int id,int s,int t,double dm,bool needpath,double OPEN);

	//dfs
	vector<int> visit;

	//独裁
	double energy;
	double delay;

public:
	CGraph(char* inputFile);

	int canNotReach(int s, int t){
		visit.resize(n, 0);
		dfs(s);
		return !visit[t];
	}

	void clearOcc(){
		for(int i=0;i<m;i++){
			Link[i]->load = 0; 
			Link[i]->use = 0;
			Link[i]->latency = INF;
		}
	}
	void SetUNVISITED(){
		for(int i = 0;i < n;++i){
			visit[i] = false;
		}
		DFSflag = 0;
		pathver.clear();
	}

	~CGraph(){
		for(unsigned int i = 0; i < Link.size(); i++)
			if(!Link[i])
				delete Link[i];
	}
};


CGraph::CGraph(char* inputFile)
{
	ifstream file(inputFile);
	file >> n >> m;
	adjL.resize(n); 
	adjRL.resize(n); 
	reqPathID.resize(200);
	set<int> vert;
	int a, b; double d;double c;
	for (int i = 0; i < m; i++)
	{
		file >> a >> b >> c >> d;
		vert.insert(a);
		vert.insert(b);
		c = rand()%MAXWEIGHT;
		CEdge *e=new CEdge(i,a,b,c,d+MINCAPACITY);
		Link.push_back(e);
		adjL[a].push_back(e); //出度边
		adjRL[b].push_back(e); //入度边
	}
	file.close();

	vector<CEdge*> emptylist;
    vector<CEdge*>::iterator it;
    for(it=Link.begin(); it!=Link.end(); it++)
    {
        mapVID_Vertex.insert(pair<int,CVertex*>((**it).getTail(),new CVertex((**it).getTail())));
        mapVID_Vertex.insert(pair<int,CVertex*>((**it).getHead(),new CVertex((**it).getHead())));
        mapVID_listEdge.insert(pair<int,vector<CEdge*>>((**it).getTail(),emptylist));
        mapVID_listEdge.insert(pair<int,vector<CEdge*>>((**it).getHead(),emptylist));
        mapVID_listEdge[(**it).getTail()].push_back(*it);
    }

	set<int>::iterator i;
	for(i=vert.begin();i!=vert.end();i++){ 
		ver.push_back(*i);
	}

}

double CGraph::dijkstra(int id,int s, int t, double dm ,bool ORselfish,bool needpath){
	vector<int> p, flag;
	vector<double> d;//带宽利用率
	for(int i = 0; i < n; i++){
		p.push_back(-1);
		flag.push_back(0);
		d.push_back(INF);
	}

	if(ORselfish){
		for(int i = 0; i < m; i++){
			CEdge *e = Link[i];
			if(e->capacity - e->use <= dm)
				e->latency = INF;
			else
				e->latency = 1.0/(e->capacity - e->use - dm);
		}
	}

	d[s] = 0;
	int cur = s;
	do{	
		flag[cur] = 1;
		for(unsigned int i = 0; i < adjL[cur].size(); i++){
			CEdge *e = adjL[cur][i];
			if(CONSTANT){  
				if(e->capacity - e->use >= dm && d[e->head] > d[e->tail] + e->dist){
					d[e->head] = d[e->tail] + e->dist;
					p[e->head] = e->id;
				}
			}
			else {
				if(ORselfish){  //latency	
					if( d[e->head] > d[e->tail] + e->latency ){ 
						d[e->head] = d[e->tail] + e->latency;
						p[e->head] = e->id;
					}		
				}
				else{  //util	
					double util = ( dm + e->use )/e->capacity; // 带宽利用率最小
					double tail_util = max(util,d[e->tail]);
					if(e->capacity - e->use >= dm && d[e->head] > tail_util ){
						d[e->head] = tail_util;
						p[e->head] = e->id;
					}
				}
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

double CGraph::dijkstraOR(int s, int t,double dm){
	vector<int> p, flag;
	vector<double> d;//latency
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
			if( d[e->head] > d[e->tail] + e->latency ){ ////latency	
				d[e->head] = d[e->tail] + e->latency;
				p[e->head] = e->id;
			}		
		}
		cur = -1;
		for(unsigned int i = 0; i < this->ver.size(); i++)
			if(!flag[ver[i]] && (cur == -1 || d[cur] > d[ver[i]] ))
				cur = ver[i];
	}while(cur != -1);

	cur = t;
	do{
		if(p[cur] == -1)
			break;
		Link[p[cur]]->use += dm;
		cur = Link[p[cur]]->tail;
	}while(cur != s);

	return d[t];
}

void genGraph(int n, int m, char route[]){ 
	FILE *out = fopen(route, "w");
	fprintf(out,"%d %d\n",n,m);
	for(int i = 1; i < min(n, m+1); i++){
		int t = rand()%i, w = rand()%MAXWEIGHT+1;
		int c = rand()%(MAXCAPACITY-MINCAPACITY)+2;
		fprintf(out, "%d %d %d %d\n", i, t, w,c);
	}
	for(int i = 0; i < m-n+1; i++){
		int s = rand()%n, t = rand()%n, w = rand()%MAXWEIGHT+1;
		int c = rand()%(MAXCAPACITY-MINCAPACITY)+2;
		while(t == s)
			t = rand()%n;
		fprintf(out, "%d %d %d %d\n", t, s, w, c);
	}
	fclose(out);
}

void genGraphOR(int n1,unsigned int n2,int m,char route[]){ 
	FILE *out=fopen(route,"w");
	fprintf(out,"%d %d\n",n1,m);
	set<int> ver; //不存重复元素
	while(ver.size()<n2){
		int s=rand()%n1;
		ver.insert(s);
	}

	vector<int> ver2;
	set<int>::iterator i=ver.begin();
	for(;i!=ver.end();i++)
		ver2.push_back(*i);

	for(unsigned int i = 1; i < n2; i++){
		int t = rand()%i, w = rand()%MAXWEIGHT+1;
		int c = rand()%(MAXCAPACITY-MINCAPACITY)+2;
		fprintf(out, "%d %d %d %d\n", ver2[i], ver2[t], w, c);
	}

	for(unsigned int i = 0; i < m-n2+1; i++){
		int s = rand()%n2, t = rand()%n2, w = rand()%MAXWEIGHT+1;
		int c = rand()%(MAXCAPACITY-MINCAPACITY)+2;
		while(t == s)
			t = rand()%n2;
		fprintf(out, "%d %d %d %d\n", ver2[s], ver2[t], w, c);
	}
	fclose(out);
}

double linearCal(double x, double c)
{
	double util=x/c;
	if(util<0.3333) return x*1;
	else if(util<0.6667) return x*3-2.0*c/3.0;
	else if(util<0.9) return x*10-16*c/3;
	else return x*70-178*c/3;
}

#endif