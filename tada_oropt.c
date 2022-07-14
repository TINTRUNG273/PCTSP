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

void show_array(int tour[MAX_N], int len) {
  int i;
  printf("[");
  for(i=0;i<len-1;i++) {
    printf("%d ",tour[i]);
  }
 printf("%d]\n",tour[i]);
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
  int count=1,countPrec=0,finish=0;
  double min_d, dist1, dist2;

  while(count != 0){
    count=0;
    min_d=0;

    for(i=0;i<=n-3;i++){
      j=i+1;
      for(k=i+2;k<=n-1;k++){
        l=(k+1)%n;
        a = tour[i];
        b = tour[j];
        c = tour[k];
        d = tour[l];
        countPrec = 0;
        finish = 0;
        for(int r=j;r<=k;r++){
    			for(int z=0;z<m;z++){
    				if(tour[r]==prec[z])
    					countPrec++;
    				if(countPrec>1){
    					finish++;
    					break;
    				}
    			}
    			if(finish>0) break;
    		}
  		  if(finish>0)
  			  continue;

        if(dist(p[a],p[b])+dist(p[c],p[d])>dist(p[a],p[c])+dist(p[b],p[d])){
          if(min_d < dist(p[a],p[b])+dist(p[c],p[d])-dist(p[a],p[c])-dist(p[b],p[d])){
            min_d = dist(p[a],p[b])+dist(p[c],p[d])-dist(p[a],p[c])-dist(p[b],p[d]);
            g=j;
            h=k;
            count+=1;
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
  printf("2-opt changed = %lf\n", tour_length(p,n,tour));
//  sprintf(tourFileName,  "tour%08d.dat",++num);
//  write_tour_data(tourFileName,n,tour);
}

void revcpy(int NewTour[], int tour[], int n, int dst, int src, int num) {
    for (int i = 0; i < num; i++) {
        NewTour[(dst+i)%n] = tour[(src-i)%n];
    }
}
/*
void ThreeOpt(struct point p[MAX_N],  int n, int tour[],int m, int prec[MAX_N]) {
    int changed = 1;
    int newTour[MAX_N];
    int a,b,c,d,e,f;
    int i, j, k, mm;
    int min;
    double dd[5];
    int countPrec=0,finish=0;

    while (changed) {
        changed = 0;
        for (i = 0; i < n-5; i++) {
            for (j = i+1; j < n-3; j++) {
                for (k = j+1; k < n-1; k++) {
                    min = 0;
                    a=tour[i]; b=tour[i+1];
                    c=tour[j]; d=tour[j+1];
                    e=tour[k]; f=tour[k+1];

                    countPrec = 0;
                    finish = 0;
                    for(int r=j;r<=k;r++){
                        for(int z=0;z<m;z++){
                            if(tour[r]==prec[z])
                                countPrec++;
                            if(countPrec>1){
                                finish++;
                                break;
                            }
                        }
                        if(finish>0) break;
                    }
                    if(finish>0)
                        continue;

                    dd[0] = dist(p[a],p[b]) + dist(p[c],p[d]) + dist(p[e],p[f]);  // 元の距離
                    dd[1] = dist(p[a],p[e]) + dist(p[d],p[b]) + dist(p[f],p[c]);
                    dd[2] = dist(p[d],p[a]) + dist(p[e],p[c]) + dist(p[f],p[b]);
                    dd[3] = dist(p[c],p[f]) + dist(p[a],p[d]) + dist(p[e],p[b]);
                    dd[4] = dist(p[a],p[c]) + dist(p[d],p[f]) + dist(p[e],p[b]);
                    
                    for (mm = 1; mm <= 4; mm++) {
                        if (dd[mm] < dd[min]) min = mm;    // ４通りのつなぎ替えを比較
                    }
                    if (min == 0) continue;

                    // 繋ぎ替えを実行
                    memcpy(&newTour[0], &tour[i+1], (j-i)*sizeof(int)); // i+1～j
                    if (min == 1) { // D(i,k) + D(j+1,i+1) + D(k+1,j)
                        memcpy(&newTour[j-i], &tour[k+1], (n-k-1)*sizeof(int)); // k+1～N
                        memcpy(&newTour[j-i+n-k-1], &tour[0], (i+1)*sizeof(int));   // 0～i
                        revcpy(newTour, tour, n, (j+n-k)%n, k, k-j);    // k～j+1
                    } else if (min == 2) { // D(j+1,i) + D(k,j) + D(k+1,i+1);
                        revcpy(newTour, tour, n, j-i, k, k-j);  // k～j+1
                        revcpy(newTour, tour, n, k-i, n+i, n+i-k);  // k+1～i
                    } else if (min == 3) { // D(j,k+1) + D(i,j+1) + D(k,i+1)
                        memcpy(&newTour[j-i], &tour[k+1], (n-k-1)*sizeof(int)); // k+1～N
                        memcpy(&newTour[j-i+n-k-1], &tour[0], (i+1)*sizeof(int));   // 0～i
                        memcpy(&newTour[j+n-k], &tour[j+1], (k-j)*sizeof(int)); // j+1～k
                    } else if (min == 4) { // D(i,j) + D(j+1,k+1) + D(k,i+1)
                        revcpy(newTour, tour, n, j-i, n+i, n+i-k);  // k+1～i
                        memcpy(&newTour[j+n-k], &tour[j+1], (k-j)*sizeof(int)); // j+1～k
                    }
                    //printf("3-opt improved! min=%d\n", min);
                    memcpy(tour, newTour, n*sizeof(int));
                    changed = 1;
                    
                }
            }
        }
        printf("%d",changed);
    }
  printf("3-opt changed = %lf\n", tour_length(p,n,tour));
	sprintf(tourFileName,  "tour%08d.dat",++num);
	write_tour_data(tourFileName,n,tour);   
}
*/

void ni(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]){
FILE *fp;
  int i,j,a,b,c,r;
  int visited[MAX_N]; // 都市iが訪問済みならば1そうでなければ0
  double d[MAX_N]; // 未訪問点 r から現在の部分巡回路までの最短距離を d[r] に保存
  int nearest[MAX_N]; /* 未訪問点 r を現在の部分巡回路内の枝(i,i+1)に挿入する
                        ときに最も距離の増加が小さい i を nearest[r]に保存*/
  double dist1,dist2, min_dist;
  int min_i,min_j,min_r;
  int sbtlen=0;
  int er, l, k = 2;

  for(a=0;a<n;a++) visited[a]=0; // 最初は全都市は未訪問

  // a= 0 に最も近い点を探す
  for(i = 0; i < m; i++) {
    tour[i]=prec[i];
    visited[prec[i]] =1;

  }
  sbtlen = m;

  while(sbtlen<n) {
    min_dist=INF;
    for(r=0;r<n;r++) {
      if(!visited[r]) {
	      d[r]=INF;
          for(l = k; l < m; l++) {
        if(prec[l] == r) er = 1;
        }
        if(er == 1) continue;
	      for(i=0;i<sbtlen;i++) {
	        a=tour[i];
          if(dist(p[r],p[a])<d[r]) {
            nearest[r]=i;
            d[r]=dist(p[r],p[a]);
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
      a=tour[i];
      if(i-1>=0){
        min_i=i-1;
        b=tour[min_i];
      } else {
        min_i=sbtlen-1;
        b=tour[min_i];
      }
      j = (i+1)%sbtlen;
      c=tour[j];
	    if (dist(p[r],p[a])+dist(p[r],p[b])-dist(p[a],p[b])<dist(p[r],p[c])+dist(p[r],p[a])-dist(p[a],p[c])) {
	      insert(tour,&sbtlen,min_i+1,r);
	    } else {
        insert(tour,&sbtlen,j,r);
      }
    
    //    printf("r,i,j,d[r] = %d %d %d %lf\n", r,i,j,d[r]);
    visited[r]=1;
    if(prec[k-1] == nearest[r]) k++;
    
  }
  printf("ni changed = %lf\n", tour_length(p,n,tour));
}

void nn(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]){ //自分たちで書く！！！！
  int i,j,r,nearest=-1,s,z;
  int a,b,count=0;
  int visited[MAX_N]; // 都市iが訪問済みならば1そうでなければ0
  int ttour[MAX_N];
  double min,mmin=INF;

  // visitedの中身
  // 1　訪問
  // -1 未訪問(制約条件内)
  // 0 未訪問

  for(int k=0;k<n;k++){

    for(i=0;i<n;i++)
      ttour[i]=0;

    for(r=0;r<n;r++)
      visited[r]=0; // 最初は全都市は未訪問
    
    for(r=1;r<m;r++)
      visited[prec[r]]=-1; //制約条件の最初の都市以外の都市の未訪問を区別する

    // for(r=0;r<n;r++)
    //   printf("%d ",visited[r]);
    
    // printf("\n");

    if(visited[k]==-1) continue;
    
    r=0; //このrが制約条件の配列の位置を表す
    if(k==prec[r]) r++;

    ttour[0]=k;         // 最初に訪問する都市は k
    visited[k]=1;      // 都市kは訪問済み
  
    for(i=0;i<n;i++) {
        min=INF;
        a = ttour[i];

        //最後に訪問した都市 a == tour[i]から最短距離にある未訪問都市nearestを見つける
        for(j=0;j<n;j++){
            if(visited[j] != 1 && visited[j] != -1){
              if(dist(p[a],p[j])<min){
                  nearest=j;
                  min=dist(p[a],p[j]);
              }
            }
        }
      
        ttour[i+1]=nearest; // i+1 番目に訪問する都市を nearest にして, 
        visited[nearest]=1;// nearest を訪問済みとする.
    
        if(nearest==prec[r]){
          r += 1;
          visited[prec[r]]=0;
        }

        count = 0;

        for(z=0;z<n;z++){ //制約条件以外の未到達点がないか探す
          if(visited[z]==0){
            count++;
          }
        }
        if(count == 0){
          i++;
          break;
        }
    }
     if(count==0){
      for(r=r;r<m;r++){
        i+=1;
        ttour[i]=prec[r];
      }
    }
    // printf("tour   :"); show_array(ttour,n);
    // printf("visited:"); show_array(visited,n);
    if(mmin>tour_length(p,n,ttour)){
        mmin=tour_length(p,n,ttour);
        for(i=0;i<n;i++) tour[i]=ttour[i];
        //printf("tour   :"); show_array(ttour,n);
        // printf("visited:"); show_array(visited,n);
    }
  }
  printf("nn changed = %lf\n", tour_length(p,n,tour));
}

void ci(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]){
  int i,j,a,b,r;
  int visited[MAX_N]; // 都市iが訪問済みならば1そうでなければ0
  double d[MAX_N]; // 未訪問点 r から現在の部分巡回路までの最短距離を d[r] に保存
  int nearest[MAX_N]; /* 未訪問点 r を現在の部分巡回路内の枝(i,i+1)に挿入する
                        ときに最も距離の増加が小さい i を nearest[r]に保存*/
  double dist1,dist2, min_dist;
  int min_i,min_j,min_r;
  int sbtlen=0;
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
    if(prec[i]==nearest[r]) i++;
   
  
  }
  printf("ci changed = %lf\n", tour_length(p,n,tour));
  
}

void orOpt1(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]){
  int a,b,c,d,e;
  int i0, i, i1, j, j1, g;
  int tmp, count=1;
  double l1, l2 ,l3, l4, l5, l6;
  
	while(count>0){
		count=0;
		for(i=0;i<n;i++){
		  if (check_prec(prec, tour[i] , m) > 0)  {
		    continue;
		  }
			i0 = i - 1;
			if(i0 < 0) i0 = n - 1;
			i1 = (i + 1) % n;
			
			for(j=i;j<n;j++){
			  if (check_prec(prec, tour[i] , m) > 0) continue;
			  j1 = (j + 1) % n;
			  a = tour[i0]; 
        b = tour[i1];
        c = tour[i]; 
        d = tour[j]; 
        e = tour[j1]; 
				if(j != i && j1 != i ){
				  l1 = dist(p[a],p[c]);
				  l2 = dist(p[c],p[b]);
				  l3 = dist(p[d],p[e]);
				  l4 = dist(p[a],p[b]);
			    l5 = dist(p[d],p[c]);
				  l6 = dist(p[c],p[e]);

				  if(l1 + l2 + l3 > l4 + l5 + l6){
					  count++;
					  tmp = tour[i];
					  for(g = i; g < j; g++) tour[g] = tour[g+1];
					  tour[j] = tmp;    				    
					}
				}
		  }		
	  }
//		sprintf(tourFileName,  "tour%08d.dat",++num);
//		write_tour_data(tourFileName,n,tour);
		printf("orOpt1 changed = %lf\n", tour_length(p,n,tour));
		
  }
}

// ここから下が大体やったところ

void orOpt2(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]){
  int a,b,c,d,e;
  int i0, i, i1, j, j1, g;
  int tmp, count=1;
  double l1, l2 ,l3, l4, l5, l6;
  
	while(count>0){
		count=0;
		for(i=0;i<n;i++){
		  if (check_prec(prec, tour[i] , m) > 0)  {
		    continue;
		  }
			i0 = i - 1;
			if(i0 < 0) i0 = n - 1;
			i1 = (i + 1) % n;
			
			for(j=i;j<n;j++){
			  if (check_prec(prec, tour[i] , m) > 0) continue;
			  j1 = (j + 1) % n;
			  a = tour[i0]; 
        b = tour[i1];
        c = tour[i]; 
        d = tour[j]; 
        e = tour[j1]; 
				if(j != i && j1 != i ){
				  l1 = dist(p[a],p[c]);
				  l2 = dist(p[c],p[b]);
				  l3 = dist(p[d],p[e]);
				  l4 = dist(p[a],p[b]);
			    l5 = dist(p[d],p[c]);
				  l6 = dist(p[c],p[e]);

				  if(l1 + l2 + l3 > l4 + l5 + l6){
					  count++;
					  tmp = tour[i];
					  for(g = i; g < j; g++) tour[g] = tour[g+1];
					  tour[j] = tmp;    				    
					}
				}
		  }		
	  }
//		sprintf(tourFileName,  "tour%08d.dat",++num);
//		write_tour_data(tourFileName,n,tour);
		printf("orOpt1 changed = %lf\n", tour_length(p,n,tour));
		
  }
}

void orOpt_k(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N], int k){
  int i, j, l;
  int a, b, c, d, e, f;
  int curi,curj;
  double tourMin=INF, curMin=INF, distMin=INF;
  double dist1,dist2;
  int changed=0, mode;
  int prec_s, prec_e, prec_num, prec_mode, prec_1, prec_2;
  int max_n=n;

/**/
  //ついでに
  int prec_tour[MAX_N];
  for(l=0;l<n;l++){
    for(j=0;j<m;j++){
      if(tour[l]==prec[j]){
        prec_tour[l]=j;
        break;
      }else{
        prec_tour[l]=-1;
      }
    }
  }


//------
/*
  prec_1=check_prec_tx(k,prec_tour,n,m,prec,1);         // ここら辺まとめられるかも
  prec_2=check_prec_tx(k,prec_tour,n,m,prec,2);
  if(prec_1 == m-1){
    if(prec_2 == 0){
      prec_mode=0;      // 通常の順番
    }else{
      prec_mode=1;      // 逆順
    }
  }else{
    if(prec_1 < prec_2){
      prec_mode=0;
    }else{
      prec_mode=1;
    }
  }
*/
/**/
//-------


//  for(i=0;i<max_n;i++){
  for(i=0;i<n;i++){               // 比較するk近傍の始点
//    prec_num=check_prec_num(k,tour,i,m,prec);         // 塊の中のprecの数
    changed=0;
    prec_num=0;
    for(j=i;j<i+k;j++){                                 // 塊の中のprecの数
      if(prec_tour[j]!=-1){
        prec_num++;
        prec_e=prec_tour[j];
      }
    }
    
    if(prec_num==0 || prec_num==m-1 || prec_num==m){  // 0,m-1,m個なら通常通り(どこでも可)
      mode=0;
    }else{                                            // それ以外なら場合分け
//      prec_s=check_prec_start(k,tour,i,m,prec);       // 塊の中のprecの最初
      prec_s=(prec_e-prec_num+1+m)%m;
//      prec_e=(prec_s+prec_num-1)%m;                     // 塊の中のprecの最後
/*      if(prec_e<prec_s){
        mode=2;
      }else{
        mode=1;
      }
*/
      mode=1;
    }
    distMin=INF;
    for(j=i+k;j+1<n+i;j++){         // 間に入れるところ(最初を比較する塊の次からにする)
      if(mode==1){                  // 入れるところが制限に引っかかるなら引っかからないところまで休止
        if(tour[j%n]==prec[(prec_e+1)%m]){
          mode=3;
          continue;
        }
//      }else if(mode==2){

      }else if(mode==3){            // 待機
        if(tour[j%n]==prec[(prec_s-1+m)%m]){
          mode=0;
        }else{
          continue;
        }
      }
      a=tour[j%n];        b=tour[(j+1)%n];          // 比較するところ
      c=tour[(i-1+n)%n];  d=tour[i];                // 塊の始点
      e=tour[(i+k-1)%n];  f=tour[(i+k)%n];          // 塊の終点
      dist1=dist(p[a],p[b])+dist(p[c],p[d])+dist(p[e],p[f]);        // 変える前
      dist2=dist(p[a],p[d])+dist(p[b],p[e])+dist(p[c],p[f]);        // 変えた後

      if(dist1>dist2){
        if(dist2-dist1<distMin){
          distMin=dist2-dist1;
          curi=i;  curj=j;
          changed=1;
        }
      }
    }

    if(changed){
//      printf("i:%d, j:%d\n", i, curj);
      change_tour(n, tour, k, curi, curj, p);
//      show_array(tour,n);
      printf("len:%lf\n",tour_length(p,n,tour));
      
//      printf("len:%lf\n",tour_length(p,n,tour));
//      max_n-=k;
      for(l=0;l<n;l++){
        for(j=0;j<m;j++){
          if(tour[l]==prec[j]){
            prec_tour[l]=j;
            break;
          }else{
            prec_tour[l]=-1;
          }
        }
      }
    }
  }
}

void change_tour(int n, int tour[MAX_N], int k, int curi, int curj, struct point p[MAX_N]){
  int i, j, s=0;
  int curTour[MAX_N];
  int mode;
  if(curi+k>=n){
    mode=1;
  }else{
    mode=0;
  }
  if(mode==0){
    for(i=0;i<n;i++){               // curTourにOrOptで変えたtourを入れる
      if(curi<=(curj+1+i)%n && (curj+1+i)%n<curi+k){
      }else{
        curTour[s]=tour[(curj+1+i)%n];
        s++;
      }
    }
    for(i=0;i<k;i++){
      curTour[n-k+i]=tour[curi+i];
    }
  }else if(mode==1){
    for(i=0;i<n;i++){               // curTourにOrOptで変えたtourを入れる
      if((curi+k)%n<=(curj+1+i)%n && (curj+1+i)%n<curi){
        curTour[s]=tour[(curj+1+i)%n];
        s++;
      }else{
      }
    }
    for(i=0;i<k;i++){
      curTour[n-k+i]=tour[curi+i];
    }
  }
  
  for(i=0;i<n;i++){
    tour[i]=curTour[i];
  }
//  show_array2(n,tour);
//  printf("len:%lf\n",tour_length(p,n,tour));
}

int check_prec_start(int k, int tour[MAX_N], int n, int i, int m, int prec[MAX_N]){
  int p, q;
  for(p=0;p<k;p++){       // 塊 の最初のprec を返す
    for(q=0;q<m;q++){     // prec にあるか あるならその場所ｑを返す
      if(tour[(i+p)%n]==prec[q]) return q;
    }
  }
  return -1;
}

int check_prec_num(int k, int tour[MAX_N], int n, int i, int m, int prec[MAX_N]){
  int p, q, num=0;
  for(p=0;p<k;p++){       // 塊 の最後のprec を返す
    for(q=0;q<m;q++){     // prec にあるか あるならその場所ｑを返す 最後まで探す
      if(tour[(i+p)%n]==prec[q]){
        num++;
        break;
      }
    }
  }
  return num;
}

int check_prec_x(int k, int tour[MAX_N], int n, int i, int m, int prec[MAX_N], int x){
  int p, q, num=0;
  for(p=0;p<k;p++){       // 塊 の最後のprec を返す
    for(q=0;q<m;q++){     // prec にあるか あるならその場所ｑを返す 最後まで探す
      if(tour[(i+p)%n]==prec[q]){
        num++;
        if(num==x) return q;
        break;
      }
    }
  }
  return -1;
}
int check_prec_tx(int k, int prec_tour[MAX_N], int n, int m, int prec[MAX_N], int x){
  int p, num=0;
  for(p=0;p<k;p++){       // x番目のprec_tourの値 を返す
    if(prec_tour[p]!=-1){
      num++;
      if(num==x) return prec_tour[p];
    }
  }
  return -1;
}

int check_prec2(int value, int m, int prec[MAX_N]){
  int i;
  for(i=0;i<m;i++){
    if(value==prec[i]){
      return i;
    }
  }
  return -1;
}

void show_array2(int n, int tour[MAX_N]){
  int i;
  for(i=0;i<n;i++){
    printf("tour[%d]:%d\n",i,tour[i]);
  }
}

int main(int argc, char *argv[]) {
  int n;                   // 点の数 
  int m;                   // 順序制約に現れる点の数
  struct point p[MAX_N];   // 各点の座標を表す配列 
  int tour[MAX_N];   // 巡回路を表現する配列
  int prec[MAX_N];   // 順序制約を表現する配列
  int i,j,k;

  if(argc != 2) {
    fprintf(stderr,"Usage: %s <tsp_filename>\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  // 点の数と各点の座標を1番目のコマンドライン引数で指定されたファイルから読み込む
  read_tsp_data(argv[1],p,&n,prec,&m);

  //順序制約の確認
  //for(i=0;i<m;i++) printf("%d\n",prec[i]);

  // 最近近傍法による巡回路構築
  ci(p,n,tour,m,prec);
//  show_array2(n,tour);
  /*
  for(i = 0; i < 100; i++){
    printf("i = %2d : \n",i);
    s = rand() % n;
        if(check_prec(prec, s , m) > 0)  
          while(check_prec(prec, s , m) > 0)  s = rand() % n;
        ci(p,n,tour,m,prec);
  }
  */
  //ci(p,n,tour,m,prec); 
//  TwoOpt(p,n,tour,m,prec);
  //ThreeOpt(p,n,tour,m,prec);
//  orOpt1(p,n,tour,m,prec); 
  //TwoOpt(p,n,tour,m,prec);
//  show_array(tour,n);
  for(k=5;k>0;k--){
    orOpt_k(p,n,tour,m,prec,k);
  }
  // ファイルに出力
  write_tour_data("tour1.dat",n,tour);
  printf("%lf\n",tour_length(p,n,tour)); // 巡回路tourの長さ
  
  exit(EXIT_SUCCESS);
}