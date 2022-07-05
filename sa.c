#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define MAX_N 10000 
#define INF 100000000
#define EPSILON 0.00000001
#define SWAP(a,b){int temp; temp=(a); (a)=(b); (b)=temp; }

struct CELL{
  int data;
  struct CELL *next;
};

struct point {
  int x;
  int y;
};

double dist(struct point p, struct point q) {  
  return sqrt((p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y));
}

double tour_length(struct point p[MAX_N], int n, int tour[MAX_N]) {
  int i;
  double sum=0.0;
  for(i=0;i<n;i++) sum+=dist(p[tour[i]],p[tour[(i+1)%n]]);
  return sum;
}

void read_tsp_data(char *filename, struct point p[MAX_N],int *np) {
  FILE *fp;
  char buff[100];
  int i;

  if ((fp=fopen(filename,"rt")) == NULL) {// 指定ファイルを読み込み用に開く
    fprintf(stderr,"Error: File %s open failed.\n",filename);
    exit(0);
  }   

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // DIMENSION で始まる行に出会う
	&&(strncmp("DIMENSION",buff,9)!=0)) ; // まで読み飛ばす. 
  sscanf(buff,"DIMENSION : %d",np);           // 点の数 *np を読み込む

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // NODE_COORD_SECTIONで始まる
	&&(strncmp("NODE_COORD_SECTION",buff,18)!=0)) ; // 行に出会うまで, 
                                                        // 読み飛ばす. 
  for(i=0;i<*np;i++) {                       // i=0 から i=(*np)-1まで
    if(fgets(buff,sizeof(buff),fp)!=NULL) 
      sscanf(buff,"%*d %d %d",&(p[i].x),&(p[i].y)); // i番目の点の座標を読み込む
  }                                 

  fclose(fp);
}

double temperature(double alpha, double s_t){
  int i;
  double t;
  
  s_t = alpha*s_t;
  
  return s_t;
}

double probability(double best_length, double n_length, double t){
  if(n_length < best_length){
    return 1;
  }else{
    return exp((best_length-n_length)/t);
  }
}

void sa(struct point p[MAX_N], int n, int tour[MAX_N]){
  double s_t=100, e_t=0.1;
  int i,j,t,k,x,ransuu1,ransuu2;
  double kakuritsu;
  double alpha = 0.99, max_step = 1000, step;

  double kaunt=0,kaunt2=0;

  int best_tour[MAX_N], n_tour[MAX_N], t_tour[MAX_N],c_tour[MAX_N],visited[MAX_N];
  double best_length, n_length, t_length, c_length;
    
    
  srand((unsigned int)time(NULL));
  
  for(i=0; i<n; i++){
    visited[i] = 0;
  }
  i = 0;
  do{
    j = rand()%n;
    if(visited[j] == 0){
      tour[i] = j;
      visited[j] = 1;
      i++;
    }
  }while(i != n);
  
  for(i=0; i<n; i++){
    best_tour[i] = tour[i];
  }
  best_length = tour_length(p,n,best_tour);
  
  for(i=0; i<n; i++){
    n_tour[i] = best_tour[i];
  }
   
  
  do{
    printf("%f  %f\n", s_t, tour_length(p,n,best_tour));
    for(step=0; step<max_step; step++){
      
      c_length = INF;
    
      for(x=0; x<10; x++){
        ransuu1 = rand()%(n-1);
        do{
          ransuu2 = rand()%(n-1);
        }while(ransuu1 == ransuu2);
    
        //2-opt
        if(ransuu1 < ransuu2){
          while(ransuu1 < ransuu2){
            SWAP(n_tour[ransuu1], n_tour[ransuu2]);
            ransuu1++; ransuu2--;
          } 
        }else{
          while(ransuu1 > ransuu2){
            SWAP(n_tour[ransuu2], n_tour[ransuu1]);
            ransuu2++; ransuu1--;
          }
        }
      
        if(c_length > tour_length(p,n,n_tour)){  
          for(i=0; i<n; i++){
            c_tour[i] = n_tour[i];
          }
          c_length = tour_length(p,n,n_tour);
        }   
      } 

      if(c_length < best_length){
        best_length = c_length;
        for(i=0; i<n; i++){
          best_tour[i] = c_tour[i];
          n_tour[i] = c_tour[i];
        }
        kaunt++;
      }  
  
      kakuritsu = probability(best_length, c_length, s_t);
      if(rand()/RAND_MAX <= kakuritsu){
        for(i=0; i<n; i++){
          n_tour[i] = c_tour[i];
        }
        kaunt2++;
      }
      
    }
    s_t = temperature(alpha,s_t);
    
  }while(s_t>e_t);

  printf("Number of times the best tour was replaced ");
  printf("%f\n",kaunt);
  printf("Number of probability transitions ");
  printf("%f\n",kaunt2);

  for(i=0; i<n; i++){
    tour[i] = best_tour[i];
  }
}

void write_tour_data(char *filename, int n, int tour[MAX_N]){
  FILE *fp; 
  int i;
 
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

void nn(struct point p[MAX_N], int n, int tour[MAX_N]){
  int i,j,r,nearest=-1;
  int a,b;
  int visited[MAX_N]; // 都市iが訪問済みならば1そうでなければ0
  int ttour[MAX_N];
  double min,mmin=INF;

  for(int k=0;k<n;k++){

    for(r=0;r<n;r++)
        visited[r]=0; // 最初は全都市は未訪問

    ttour[0]=k;         // 最初に訪問する都市は k
    visited[k]=1;      // 都市kは訪問済み

    for(i=0;i<n;i++) {
        min=INF;
        a = ttour[i];

        //最後に訪問した都市 a == tour[i]から最短距離にある未訪問都市nearestを見つける
        for(j=0;j<n;j++){
            if(visited[j] != 1){
              if(dist(p[a],p[j])<min){
                  nearest=j;
                  min=dist(p[a],p[j]);
              }
            }
        }

        ttour[i+1]=nearest; // i+1 番目に訪問する都市を nearest にして, 
        visited[nearest]=1;// nearest を訪問済みとする. 
        // printf("tour   :"); show_array(ttour,n);
        // printf("visited:"); show_array(visited,n);

    }

    if(mmin>tour_length(p,n,ttour)){
        mmin=tour_length(p,n,ttour);
        for(i=0;i<n;i++) tour[i]=ttour[i];
        //printf("tour   :"); show_array(ttour,n);
        //printf("visited:"); show_array(visited,n);
    }

}
}

void TwoOpt(struct point p[MAX_N], int n, int tour[MAX_N]){
  int a,b,c,d;
  int i,j,k,l,g,h;
  int count=1,min_d;

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
}
void revcpy(int NewTour[], int tour[], int n, int dst, int src, int num) {
    for (int i = 0; i < num; i++) {
        NewTour[(dst+i)%n] = tour[(src-i)%n];
    }
}
void ThreeOpt(struct point p[MAX_N],  int n, int tour[]) {
    int fImproved = 1;
    int NewTour[MAX_N];
    int a,b,c,d,e,f;
    while (fImproved) {
        fImproved = 0;
        for (int i = 0; i < n-5; i++) {
            for (int j = i+1; j < n-3; j++) {
                for (int k = j+1; k < n-1; k++) {
                    int dd[5];
                    int min = 0;

                    a=tour[i];
                    b=tour[i+1];
                    c=tour[j];
                    d=tour[j+1];
                    e=tour[k];
                    f=tour[k+1];

                    dd[0] = dist(p[a],p[b]) + dist(p[c],p[d]) + dist(p[e],p[f]);  // 元の距離
                    dd[1] = dist(p[a],p[e]) + dist(p[d],p[b]) + dist(p[f],p[c]);
                    dd[2] = dist(p[d],p[a]) + dist(p[e],p[c]) + dist(p[f],p[b]);
                    dd[3] = dist(p[c],p[f]) + dist(p[a],p[d]) + dist(p[e],p[b]);
                    dd[4] = dist(p[a],p[c]) + dist(p[d],p[f]) + dist(p[e],p[b]);
                    for (int m = 1; m <= 4; m++) {
                        if (dd[m] < dd[min]) min = m;    // ４通りのつなぎ替えを比較
                    }
                    if (min == 0) continue;


                    // 繋ぎ替えを実行
                    memcpy(&NewTour[0], &tour[i+1], (j-i)*sizeof(int)); // i+1～j
                    if (min == 1) { // D(i,k) + D(j+1,i+1) + D(k+1,j)
                        memcpy(&NewTour[j-i], &tour[k+1], (n-k-1)*sizeof(int)); // k+1～N
                        memcpy(&NewTour[j-i+n-k-1], &tour[0], (i+1)*sizeof(int));   // 0～i
                        revcpy(NewTour, tour, n, (j+n-k)%n, k, k-j);    // k～j+1
                    } else if (min == 2) { // D(j+1,i) + D(k,j) + D(k+1,i+1);
                        revcpy(NewTour, tour, n, j-i, k, k-j);  // k～j+1
                        revcpy(NewTour, tour, n, k-i, n+i, n+i-k);  // k+1～i
                    } else if (min == 3) { // D(j,k+1) + D(i,j+1) + D(k,i+1)
                        memcpy(&NewTour[j-i], &tour[k+1], (n-k-1)*sizeof(int)); // k+1～N
                        memcpy(&NewTour[j-i+n-k-1], &tour[0], (i+1)*sizeof(int));   // 0～i
                        memcpy(&NewTour[j+n-k], &tour[j+1], (k-j)*sizeof(int)); // j+1～k
                    } else if (min == 4) { // D(i,j) + D(j+1,k+1) + D(k,i+1)
                        revcpy(NewTour, tour, n, j-i, n+i, n+i-k);  // k+1～i
                        memcpy(&NewTour[j+n-k], &tour[j+1], (k-j)*sizeof(int)); // j+1～k
                    }
                    //printf("3-opt improved! min=%d\n", min);
                    memcpy(tour, NewTour, n*sizeof(int));
                    fImproved = 1;
                }
            }
        }
    }
}
int main(int argc, char *argv[]) {
  int  n;                  
  struct point  p[MAX_N];   
  int tour[MAX_N];   

  if(argc != 2) {
    fprintf(stderr,"Usage: %s <tsp_filename>\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  read_tsp_data(argv[1],p, &n);
  //Simulated Annealing
  
  sa(p,n,tour);
  
  printf("%lf\n",tour_length(p,n,tour));
  nn(p,n,tour);
  TwoOpt(p,n,tour);
  ThreeOpt(p,n,tour);
  printf("%lf\n",tour_length(p,n,tour));
  
  write_tour_data("tour.dat",n,tour);

  exit(EXIT_SUCCESS);
}