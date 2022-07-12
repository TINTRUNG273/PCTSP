// GROUP 08: MEGANEs
// プログラムコンテストの第2ステージ

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define MAX_N 10000   // 点の数の最大値
#define INF 100000000 // 無限大の定義
#define EPSILON 0.00000001 //ε 小さい正の値
#define SWAP(a,b){int temp; temp=(a); (a)=(b); (b)=temp; }   

int num = 0;
char tourFileName[20];
int s;
double S_tour=INF;
int change=INF;

struct point {
  int x;
  int y;
};

double dist(struct point p, struct point q) { // pとq の間の距離を計算 
  return sqrt((p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y));
}

double tour_length(struct point p[MAX_N], int n, int tour[MAX_N]) {
  int i;
  double sum=0.0;
  for(i=0;i<n;i++) sum+=dist(p[tour[i]],p[tour[(i+1)%n]]);
  return sum;// 総距離が関数の戻り値
}

void read_tsp_data(char *filename, struct point p[MAX_N],int *np, int prec[MAX_N], int *mp) {
  FILE *fp;
  char buff[500];
  int i;

  if ((fp=fopen(filename,"rt")) == NULL) {// 指定ファイルを読み込み用に開く
    fprintf(stderr,"Error: File %s open failed.\n",filename);
    exit(EXIT_FAILURE);
  }   

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // PRECEDENCE_CONSTRAINTS:で始まる行に出会う
	&&(strncmp("PRECEDENCE_CONSTRAINTS:",buff,23)!=0)) ; // まで読み飛ばす. 
  if(strncmp("PRECEDENCE_CONSTRAINTS:",buff,23)==0)  {
    sscanf(buff+24,"%d",mp);
    for(i=0;i<*mp;i++) fscanf(fp,"%d ",&prec[i]);
  } else {
    fprintf(stderr,"Error: There is no precedence constraint in file %s.\n",filename);
    exit(EXIT_FAILURE);
  }

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // DIMENSION で始まる行に出会う
	&&(strncmp("DIMENSION",buff,9)!=0)) ; // まで読み飛ばす. 
  sscanf(buff,"DIMENSION: %d",np);           // 点の数 *np を読み込む

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // NODE_COORD_SECTIONで始まる
	&&(strncmp("NODE_COORD_SECTION",buff,18)!=0)) ; // 行に出会うまで, 
                                                        // 読み飛ばす. 
  for(i=0;i<*np;i++) {                       // i=0 から i=(*np)-1まで
    if(fgets(buff,sizeof(buff),fp)!=NULL) 
      sscanf(buff,"%*d %d %d",&(p[i].x),&(p[i].y)); // i番目の点の座標を読み込む
  }                                 

  fclose(fp);
}

void show_array(int a[MAX_N], int len) {
  int i;
  printf("[ ");
  for(i=0;i<len-1;i++) {
    printf("%3d ",a[i]);
  }
 printf("%3d ]\n",a[i]);
}

void write_tour_data(char *filename, int n, int tour[MAX_N]){
  FILE *fp; 
  int i;
 
 // 構築した巡回路をfilenameという名前のファイルに書き出すためにopen
  if((fp=fopen(filename,"wt"))==NULL){ 
    fprintf(stderr,"Error: File %s open failed.\n",filename);
    exit(EXIT_FAILURE);
  }
  fprintf(fp,"%d\n",n);
  for(i=0;i<n; i++){
   fprintf(fp,"%d ",tour[i]);
  }
  fprintf(fp,"\n");
  fclose(fp);
}

void insert(int tour[MAX_N], int* len, int k, int value) {
  int i;
  
  if(k<0 || k > *len) {
    printf("Error in insert: out of range\n");
  }
  for(i=*len;i>k;i--) {
    tour[i]=tour[i-1];
  }
  tour[k]=value;
  (*len)++;
}

int check_prec(int prec[MAX_N], int r, int m){
  int i;
  for (i =0; i < m; i++) if(prec[i] == r) {
    return 1;
  }
  return 0;
}

void TwoOpt(struct point p[MAX_N], int n, int tour[MAX_N], int m, int prec[MAX_N]){
  int a,b,c,d;
  int i,j,k,l,g,h,r,z;
  int countPrec=0,finish=0;
  double max_dist=0;

  while(max_dist != (-1)*INF ){
    max_dist=(-1)*INF;

    for(i=0;i<=n-3;i++){
      j=i+1;
      for(k=i+2;k<=n-1;k++){
        l=(k+1)%n;
        a = tour[i]; b = tour[j];
        c = tour[k]; d = tour[l];
        countPrec = -1;
        finish = 0;

        if(dist(p[a],p[b])+dist(p[c],p[d])>dist(p[a],p[c])+dist(p[b],p[d])
            && dist(p[a],p[b])+dist(p[c],p[d]) - (dist(p[a],p[c])+dist(p[b],p[d]))> max_dist){
          for(int r=j;r<=k && finish == 0;r++){
			      for(int z=0;z<m;z++){
                if(tour[r]==prec[z] && countPrec <= 0 ){
                        //countPrec++;
                        if(++countPrec == 1)
                            finish = 1;
                    }
                }
            }
					
        if(finish==0){
            g=j;
            h=k;
            max_dist = dist(p[a],p[b])+dist(p[c],p[d])-(dist(p[a],p[c])+dist(p[b],p[d]));
          }
        }
      }
    }
    while(g<h){
      SWAP(tour[g],tour[h]);
      g++;
      h--;
    }
  }
  if( S_tour > tour_length(p,n,tour)){
    S_tour = tour_length(p,n,tour);
    printf("Changed by 2-Opt  = %lf\n", tour_length(p,n,tour));
    sprintf(tourFileName,  "tour%08d.dat",++num);
    write_tour_data(tourFileName,n,tour);
    change++;
  }
}


void ci(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]){
  int i,j,a,b,r;
  int visited[MAX_N]; // 都市iが訪問済みならば1そうでなければ0
  double d[MAX_N]; // 未訪問点 r から現在の部分巡回路までの最短距離を d[r] に保存
  int nearest[MAX_N]; /* 未訪問点 r を現在の部分巡回路内の枝(i,i+1)に挿入する
                        ときに最も距離の増加が小さい i を nearest[r]に保存*/
  double dist1,dist2, min_dist;
  int min_i,min_j,min_r;
  int sbtlen=0,z=2;
  for(a=0;a<n;a++) visited[a]=0; // 最初は全都市は未訪問

  for(i = 0; i < m; i++) {
    tour[i]=prec[i];
    visited[prec[i]]=1;
  }

  sbtlen=m;
  
  while(sbtlen<n) {
    min_dist=INF;
    for(r=0;r<n;r++) {
      if(!visited[r] && (check_prec(prec, r , m) == 0)) {
	      d[r]=INF;
	      for(i=0;i<sbtlen;i++) {
	        a=tour[i];
	        j=(i+1)%sbtlen; 
	        b=tour[j];
	        if (dist(p[a],p[r])+dist(p[b],p[r])-dist(p[a],p[b])<d[r]) {
	          nearest[r]=i;
	          d[r]=dist(p[a],p[r])+dist(p[b],p[r])-dist(p[a],p[b]);
	        }
	      }
	      if (d[r]<min_dist) {
	        min_dist = d[r];
	        min_r=r;
	      }
      }
    }
    r=min_r;
    i=nearest[r];
    j=(i+1)%sbtlen;

    insert(tour,&sbtlen,j,r); 
    visited[r]=1;
    if(prec[z-1]==nearest[r]) z++;
  
  }
  printf("Changed by Ci = %lf\n", tour_length(p,n,tour));
  
}

void OrOpt1(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]){
  int a,b,c,d,e;
  int i0, i, i1, j, j1, g;
  int tmp, count=1;
  double k1, k2 ,k3, k4, k5, k6;
  
	while(count>0){
		count=0;
		for(i=0;i<n;i++){
          count=0;
          i0 = i - 1;
		  if(i0 < 0) i0 = n - 1;
		  i1 = (i + 1) % n;
          if (check_prec(prec, tour[i] , m) > 0) continue;
			for(j=0;j<n;j++){
			  if (check_prec(prec, tour[i] , m) > 0) continue;
			    j1 = (j + 1) % n;
			    a = tour[i0]; b = tour[i1];
                c = tour[i];  d = tour[j]; e = tour[j1]; 
				if(j != i && j1 != i ){
				  k1 = dist(p[a],p[c]); k2 = dist(p[c],p[b]);
				  k3 = dist(p[d],p[e]); k4 = dist(p[a],p[b]);
			      k5 = dist(p[d],p[c]); k6 = dist(p[c],p[e]);
				  if(k1 + k2 + k3 > k4 + k5 + k6){
					  count++;
					  tmp = tour[i];
                      if(i < j){
					    for(g = i; g < j; g++) tour[g] = tour[g+1];
					    tour[j] = tmp;
                        i--;
                      }else{
                        for(g = i; g > j1; g--) tour[g%n] = tour[(g-1)%n];
					    tour[j1] = tmp;
                      } 
                      break;   				        
					}
				}
		  }		
	  }
      if( S_tour > tour_length(p,n,tour)){
        S_tour = tour_length(p,n,tour);
		sprintf(tourFileName,  "tour%08d.dat",++num);
		write_tour_data(tourFileName,n,tour);
		printf("Changed by orOpt1 = %lf\n", tour_length(p,n,tour));
        change++;
      }
  }
}

void OrOpt2(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]){
  int a,b,c,d,e,f;
  int i0, i, i1, i2, j, j1, g, tmp, tmp1;
  int count=1;
  double k1, k2 ,k3, k4, k5, k6, k7 , k8, k9, k10;
  
	while(count>0){
		count=0;
		for(i=0;i<n;i++){
		  if (check_prec(prec, tour[i] , m) > 0 || check_prec(prec, tour[i+1] , m) > 0)  {
		    continue;
		  }
			i0 = i - 1;
			if(i0 < 0) i0 = n - 1;
			i1 = (i + 1) % n;
			i2 = (i1 +1) % n;
			
			for(j=i+2;j<n;j++){
			  if (check_prec(prec, tour[i] , m) > 0 || check_prec(prec, tour[i+1] , m) > 0) continue;
			  j1 = (j + 1) % n;
			  a = tour[i0]; b = tour[i2];
              c = tour[i]; d = tour[i1]; 
              e = tour[j]; f = tour[j1]; 
			  if(j != i && j1 != i ){
				  k1 = dist(p[a],p[c]); k2 = dist(p[d],p[b]);
                  k3 = dist(p[e],p[f]); k4 = dist(p[a],p[b]);
			      k5 = dist(p[e],p[c]); k6 = dist(p[d],p[f]);
				  k7 = dist(p[e],p[d]); k8 = dist(p[c],p[f]);
				  k9 = dist(p[a],p[d]); k10 = dist(p[c],p[b]);

				  if(k1 + k2 + k3 > k4 + k5 + k6 ){
					  count++;
					  tmp = tour[i];
					  tmp1 = tour[(i+1)%n];
					  for(g = i; g < j-1; g++) tour[g] = tour[(g+2)%n];
					  if(k5 + k6 > k7 + k8){
					    tour[j-1] = tmp1;    			
					    tour[j] = tmp;	
					  } 
					  else {
					    tour[j-1] = tmp;    			
					    tour[j] = tmp1;	 
					  }   
					}
				}
		  }		
	  }
      if( S_tour > tour_length(p,n,tour)){
        S_tour = tour_length(p,n,tour);      
		sprintf(tourFileName,  "tour%08d.dat",++num);
		write_tour_data(tourFileName,n,tour);
		printf("Changed by orOpt2 = %lf\n", tour_length(p,n,tour));
        change++;
      }
  }
}

void OrOpt3(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]){
  int a,b,c,d,e,f,g;
  int i0, i, i1, i2, i3, j, j1, tmp, tmp1, tmp2, count=1;
  double k1, k2 ,k3, k4, k5, k6, k7 , k8, k9, k10, k11, k12;
  double min = INF;
  
	while(count>0){
		count=0;
		for(i=0;i<n;i++){
		  if (check_prec(prec, tour[i] , m) > 0 || check_prec(prec, tour[i+1] , m) > 0)  {
		    continue;
		  }
			i0 = i - 1;
			if(i0 < 0) i0 = n - 1;
			i1 = (i + 1) % n;
			i2 = (i1 +1) % n;
			i3 = (i1 +2) % n;
			
			for(j=i+3;j<n;j++){
			  if ( ( check_prec(prec, tour[i] , m) > 0 || check_prec(prec, tour[i+1] , m) > 0 ) || check_prec(prec, tour[i+2] , m) > 0) continue;
			  j1 = (j + 1) % n;
			  a = tour[i0]; b = tour[i3]; 
              c = tour[i]; d = tour[i1]; 
              e = tour[i2]; f = tour[j]; g = tour[j1]; 
				if(j != i && j1 != i ){
				  k1 = dist(p[a],p[c]); k2 = dist(p[e],p[b]);
				  k3 = dist(p[f],p[g]); k4 = dist(p[a],p[b]);
			      k5 = dist(p[f],p[c]); k6 = dist(p[e],p[g]);
				  k7 = dist(p[f],p[c]); k8 = dist(p[f],p[d]);
				  k9 = dist(p[f],p[e]); k10 = dist(p[g],p[e]);
				  k11 = dist(p[g],p[d]); k12 = dist(p[g],p[e]);

				  if(k1 + k2 + k3 > k4 + k5 + k6 ){
					  count++;
					  tmp = tour[i];
					  tmp1 = tour[(i+1)%n];
					  tmp2 = tour[(i+2)%n];
					  for(g = i; g < j-2; g++) tour[g] = tour[(g+3)%n];
					  min = INF;	  
					  
					  if(min > k7 + k12){
					    min = k7 + k12;
						tour[j-2] = tmp; 
					    tour[j-1] = tmp1;    			
					    tour[j] = tmp2;
					  }
					  
					  if(min > k8 + k12){
					    min = k8 + k12;
						tour[j-2] = tmp1; 
					    tour[j-1] = tmp;    			
					    tour[j] = tmp2;	
					  }
					  if(min > k7 + k11){
					    min = k7 + k11;
						tour[j-2] = tmp; 
					    tour[j-1] = tmp2;    			
					    tour[j] = tmp1;
					  }
					  if(min > k8 + k10){
					    min = k8 + k10;
						tour[j-2] = tmp1; 
					    tour[j-1] = tmp2;    			
					    tour[j] = tmp;
					  }
					  if(min > k9 + k11){
					    min = k9 + k11;
						tour[j-2] = tmp2; 
					    tour[j-1] = tmp;    			
					    tour[j] = tmp1;	
					  }
					  if(min > k9 + k10){
					    min = k9 + k10;
						tour[j-2] = tmp2; 
					    tour[j-1] = tmp1;    			
					    tour[j] = tmp;
					  }
					}
				}
		  }		
	  }
	  if(S_tour > tour_length(p,n,tour)) {
	    change++;
		S_tour = tour_length(p,n,tour);
		sprintf(tourFileName,  "tour%08d.dat",++num);
		write_tour_data(tourFileName,n,tour);
		printf("Changed by orOpt3 = %lf\n", tour_length(p,n,tour));
	  }
  }
}



int main(int argc, char *argv[]) {
  int n;                   // 点の数 
  int m;                   // 順序制約に現れる点の数
  struct point p[MAX_N];   // 各点の座標を表す配列 
  int tour[MAX_N];   // 巡回路を表現する配列
  int prec[MAX_N];   // 順序制約を表現する配列
  int i,j;

  if(argc != 2) {
    fprintf(stderr,"Usage: %s <tsp_filename>\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  // 点の数と各点の座標を1番目のコマンドライン引数で指定されたファイルから読み込む
  read_tsp_data(argv[1],p,&n,prec,&m);

  // 最近近傍法による巡回路構築
    ci(p,n,tour,m,prec); 
    change=INF;
    //for(int l = 0; l <=3; l++){
    while(change!=0){
      change=0;
      TwoOpt(p,n,tour,m,prec);
      OrOpt3(p,n,tour,m,prec); 
      OrOpt2(p,n,tour,m,prec);
      OrOpt1(p,n,tour,m,prec);
    }
	
  // ファイルに出力
  write_tour_data("tour1.dat",n,tour);
  printf("%lf\n",tour_length(p,n,tour)); // 巡回路tourの長さ
  
  exit(EXIT_SUCCESS);
}