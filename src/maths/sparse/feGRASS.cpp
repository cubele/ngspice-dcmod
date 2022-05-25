# include "feGRASS.hpp"
using namespace std;

struct Edge{
    int node;
    int ind;
};

class Graph{
    public:
		int n;int m;
		Edge* edgelists;
		int* ptr;
        Graph(int n,int m):n(n),m(m){
        }
        Graph(){m=0;}
		void destroy(){
			free(edgelists);free(ptr);
		}
};

class DisjSet{
  public:
    int* parent;
    int* rank;
    DisjSet(int max_size){
		parent=(int*) malloc(sizeof(int)*max_size);
		rank=(int*) malloc(sizeof(int)*max_size);
		for(int i=0;i<max_size;i++) parent[i]=i;
    }
    int find(int x){
        int p=parent[x];
        if(parent[p]!=p){
        	p=find(p);
        	parent[x]=p;
		}
		return p;
    }
    void to_union(int x1, int x2){
        int f1 = find(x1);
        int f2 = find(x2);
        if(f1==f2) return;
        if (rank[f1] < rank[f2])
            parent[f1] = f2;
        else{
            parent[f2] = f1;
            if (rank[f1] == rank[f2])
                ++rank[f1];
		}
    }
    bool is_same(int e1, int e2){
        return find(e1) == find(e2);
    }
	void destroy(){
		free(parent);free(rank);
	}
};

struct dfsnode {int node, from;};
struct bfsnode {int node, layer;};

int *ai,*aj;
double* av;
int SIZE;int M;
int* edgeperm;
int *deg,*layer,*LCA;
int* insparsifier;
Graph gra;
int* fa,*vis,*ans,*parent;
double* dist,*score,*effwt;
double totstretch,avstretch;
bfsnode* bfsque;
double scale=-1.0;
int maxnode=1;
int totedges=0,treeedges;

double random(double range){
	long long bigint=(long long)(rand())*(long long)(RAND_MAX)+(long long)(rand());
	long long r=(long long) range;
	long long res=bigint%r;
	return (double) res;
}

inline int getfa(int x) { return fa[x] == x ? x : fa[x] = getfa(fa[x]); }

bool compareByEffWt(int i,int j){
	return effwt[i]>effwt[j];
}

bool compareByWt(Edge i,Edge j){
	return av[i.ind]>av[j.ind];
}

bool compareByDeg(int i,int j){
	return deg[i]>deg[j];
}

void BFS1(){
	int head=0,tail=0,node,curlayer,i=maxnode,comp=0;
	bfsque[tail]=(bfsnode){i,0};tail++;
	layer[i]=0;
	while(head<tail){
		node=bfsque[head].node;curlayer=bfsque[head].layer;head++;
		for(int k=gra.ptr[node];k<gra.ptr[node+1];k++){
			if(layer[gra.edgelists[k].node]<0){
				bfsque[tail++]=(bfsnode){gra.edgelists[k].node,curlayer+1};
				layer[gra.edgelists[k].node]=curlayer+1;
			}
		}
	}
	for(int i=0;i<SIZE;i++){
		if(layer[i]>=0) continue;
		head=0;tail=0;
		bfsque[tail]=(bfsnode){i,0};tail++;
		layer[i]=0;
		while(head<tail){
			node=bfsque[head].node;curlayer=bfsque[head].layer;head++;
			for(int k=gra.ptr[node];k<gra.ptr[node+1];k++){
				if(layer[gra.edgelists[k].node]<0){
					bfsque[tail++]=(bfsnode){gra.edgelists[k].node,curlayer+1};
					layer[gra.edgelists[k].node]=curlayer+1;
				}
			}
		}
	}
}

void ModKruscal3(){
	double maxdeg,sumlayer=1,scale=1;int i,j;
	BFS1();

	for(int l=0;l<M;l++){
		i=ai[l];j=aj[l];
		maxdeg=(double)(deg[i]>deg[j]?deg[i]:deg[j]);
		sumlayer=(double)(layer[i]+layer[j]);
		scale=log(maxdeg)/(sumlayer);
		effwt[l]=av[l]*scale;
	}
	int ptr1=1,ptr2=M;
	clock_t start1=clock();
    sort(edgeperm,edgeperm+M,compareByEffWt);
	clock_t end1=clock();

	DisjSet dis(SIZE);
	int ind; 
	int count=0;
	for(int l=0;l<M;l++){
		if(count>=SIZE-1) {break;}
		ind=edgeperm[l];
		i=ai[ind];j=aj[ind];
        if(!dis.is_same(i,j)){
            insparsifier[ind]=true;
			count++;
            dis.to_union(i,j);
        } 
    }
	treeedges=count;totedges=treeedges;
    cout<<"Kruscal end, "<<count<<" edges in spanning tree, "<<SIZE-count<<" connected components."<<endl;
	dis.destroy();
}

void BFSmark(int* mark,int i,int layer,int* visited,int rank,int nodenum){
	int head=0,tail=0,node,flag=0,curlayer;
	bfsque[tail]=(bfsnode){i,0};tail++;
	while(head<tail && tail<nodenum){
		node=bfsque[head].node;curlayer=bfsque[head].layer;head++;
		if(curlayer>=layer){break;}
		for(int k=gra.ptr[node];k<gra.ptr[node+1];k++){
			if(!insparsifier[gra.edgelists[k].ind]||visited[gra.edgelists[k].node]>=rank) continue;
			visited[gra.edgelists[k].node]=rank;
			bfsque[tail]=(bfsnode){gra.edgelists[k].node,curlayer+1};
			tail++;
		}
	}
	for(int j=0;j<tail;j++) {
		mark[bfsque[j].node]++;
	}
}

void addedges(double fac,int sigma){
	int ptr1=0,ptr2=M-1;
	for(int l=0;l<M;l++){
		if(insparsifier[l]){
			edgeperm[ptr2]=l;ptr2--;
		}
		else{
			edgeperm[ptr1]=l;ptr1++;
		}
	}

	clock_t start1=clock();
	sort(edgeperm,edgeperm+(ptr1),compareByEffWt);
	clock_t end1=clock();
	free(effwt);

	int uplink=1;
	int num=int(fac*double(SIZE));
	int nodenum=(int)(0.6/fac);
	int count=0;
	int ptr=0;int ind;int bdry=M/2;
	double wt,R;
	// mark->deg; visited->layer
	memset(deg,0,sizeof(int)*(SIZE));
	memset(layer,0,sizeof(int)*(SIZE));
	int i,j,m1,m2,m3,m4,m5,flag=1;
	while(1){
		if(count>=num) break;
		if(ptr>=bdry){
			ptr=0;uplink++;
			cout<<count<<" edges recovered"<<endl;continue;
		}
		ind=edgeperm[ptr];
		i=ai[ind];j=aj[ind];
		if(!insparsifier[ind] && (deg[i]<uplink && deg[j]<uplink) ){
		    count++;
			BFSmark(deg,i,sigma,layer,count,nodenum);
			BFSmark(deg,j,sigma,layer,count,nodenum);
			insparsifier[ind]=true;
	    }
		ptr++;
	}
    totedges+=count;
	cout<<count<<" edges recovered"<<endl;
    return;
}

void dfs(int start,Graph& Tree){
	stack<dfsnode> dfsstack;
	dfsstack.push((dfsnode){start,-1});
	int ptr,count=0;
	memset(deg,0,sizeof(int)*(SIZE));
	int node,child,from;int ind;double len;
	bool hasson;
	while(!dfsstack.empty()){
		node=dfsstack.top().node;
		from=dfsstack.top().from;
		hasson=false;
		int length=Tree.ptr[node+1]-Tree.ptr[node];
		int start=Tree.ptr[node];
		for(;deg[node]<length;){
			ptr=deg[node];
			child=Tree.edgelists[start+ptr].node;
			ind=Tree.edgelists[start+ptr].ind;
			if(deg[child]!=0 || !insparsifier[ind]){
				deg[node]++;
				continue;
			}
			hasson=true;
			len=1.0/(av[ind]);
			dist[child]=dist[node]+len;
			deg[node]++;
			dfsstack.push((dfsnode){child,node});
			break;
		}
		if(!hasson){
			dfsstack.pop();
			for(int k=gra.ptr[node];k<gra.ptr[node+1];k++){
				if(layer[gra.edgelists[k].node]){ 
					int ind=gra.edgelists[k].ind;
					LCA[ind]=getfa(gra.edgelists[k].node);
				}
			}
			fa[node]=from;
			layer[node]=1;
		}
	}
}

void tarjan(int root,Graph& Tree){
	totstretch=0;
	for(int i=0;i<SIZE;i++){
		fa[i]=i;
	}
	memset(layer,0,sizeof(int)*(SIZE));
	for(int i=0;i<SIZE;i++){
		if(layer[i]==0){
			dist[i]=0;
			dfs(i,Tree);
		}
	}
	int i,j;double scale,w,R,st,maxlayer,longest=0,score0;
	totstretch=0;
	// ans->edgeperm; stretch->effwt
	for(int l=0;l<M;l++){
		i=ai[l];
		j=aj[l];
		w=av[l];
		R=(dist[i]+dist[j]-2*dist[LCA[l]]);
		effwt[l]=w*R;
		totstretch+=effwt[l];	
	}
	avstretch=totstretch/M;
	cout<<"ave stretch = "<<avstretch<<endl;
	free(fa);
}

void calculateresistance(){
	tarjan(0,gra);
}

void initialize(){
    deg=(int*) malloc(sizeof(int)*(SIZE));memset(deg,0,sizeof(int)*(SIZE));
	layer=(int*) malloc(sizeof(int)*(SIZE));memset(layer,-1,sizeof(int)*(SIZE));
	dist=(double*) malloc(sizeof(double)*(SIZE));memset(dist,0,sizeof(double)*(SIZE));
	bfsque=(bfsnode*) malloc(sizeof(bfsnode)*(SIZE));
	fa=(int*) malloc(sizeof(int)*(SIZE));memset(fa,0,sizeof(int)*(SIZE));

	edgeperm=(int*) malloc(sizeof(int)*(M));LCA=(int*) malloc(sizeof(int)*(M));
	effwt=(double*) malloc(sizeof(double)*(M));
	memset(insparsifier,0,sizeof(int)*(M));
	for(int l=0;l<M;l++) {edgeperm[l]=l;}
}

void constructgraph(){
    gra=Graph(SIZE,M);
	int i,j;
	for(int l=0;l<M;l++){
		i=ai[l];j=aj[l];
		deg[i]++;deg[j]++;
	}
	gra.edgelists=(Edge*) malloc(sizeof(Edge)*(2*M));
	gra.ptr=(int*) malloc(sizeof(int)*(SIZE+1));
	int* ptrs=(int*) malloc(sizeof(int)*(SIZE));
	maxnode=0;gra.ptr[0]=0;ptrs[0]=0;
	int sum=0;
	for(int l=1;l<SIZE;l++){
		if(deg[l]>deg[maxnode]) maxnode=l;
		sum+=deg[l-1];
		gra.ptr[l]=sum;
		ptrs[l]=sum;
	}
	sum+=deg[SIZE-1];
	gra.ptr[SIZE]=sum;
	for(int l=0;l<M;l++){
		i=ai[l];j=aj[l];
		gra.edgelists[ptrs[i]]=(Edge){j,l};
		gra.edgelists[ptrs[j]]=(Edge){i,l};
		ptrs[i]++;ptrs[j]++;
	}
	free(ptrs);
}

void freeMemory(){
	gra.destroy();
	free(edgeperm);free(LCA);
	free(deg);free(layer);free(bfsque);free(dist);
}

void feGRASS(int* ai_in,int* aj_in,double* av_in,int M_in,int N_in,int* insparsifier_in,double alpha_in){
	ai=ai_in;aj=aj_in;av=av_in;M=M_in;SIZE=N_in;insparsifier=insparsifier_in;double fac=alpha_in;
	initialize();
    constructgraph();
	std::cout<<SIZE<<" nodes, "<<M<<" edges."<<std::endl;
    std::cout<<"Graph Construction finished. Graph sparsification begin..."<<std::endl;
	clock_t start=clock();
	ModKruscal3();
	calculateresistance();
	addedges(fac,1);
	clock_t end=clock();
	cout<<"Graph sparsification end. Sparsifier construction time = "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
	freeMemory();
	return;
}