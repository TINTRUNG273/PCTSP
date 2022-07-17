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
int change=INF, max_i = INF, flags[MAX_N];
double min_distance = INF;

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

// delete element at position k in array
// 配列の位置kにある要素を削除する
void delete (int tour[], int *len, int k) {
    int i;
    if (k < 0 || k > *len) {
        printf("Error in delete: out of range\n");
    }

    for (i = k; i < *len; ++i) {
        tour[i] = tour[i + 1];
    }
    (*len)--;
}

void initialize(int prec[], int m) {
    for (int i = 0; i < m; ++i) {
        int id_node = prec[i];
        flags[id_node] = 1;
    }
    sprintf(tourFileName, "tour%08d.dat", 1);
}


int check_prec(int prec[MAX_N], int r, int m){
  //int i;
  //for (i =0; i < m; i++) if(prec[i] == r) {
  //  return 1;
  //}
  //return 0;
  return flags[r];
}

// ランダム関数
//参考文献：http://c-lang.sevendays-study.com/ex-day1_macro.html
int rand_comparison(const void *a, const void *b) {
    (void)a;
    (void)b;
    return rand() % 2 ? +1 : -1;
}

// クイックソート関数でシャッフル
// 参考文献：https://monozukuri-c.com/langc-funclist-qsort/
void shuffle(void *base, size_t nmemb, size_t size) {
    qsort(base, nmemb, size, rand_comparison);
}

// get substring from string s from position start to position end
// 位置開始から位置終了までの文字列sから部分文字列を取得する
char *substr(const char *s, int start, int end) {
    char *subtext = malloc(100 * sizeof(char));
    if (start < 0)
        start = 0;
    if (end < 0) {
        end += strlen(s);
    }
    if (end < 0)
        end = 0;
    strncpy(subtext, s + start, end - start);
    subtext[end - start] = '\0';
    return subtext;
}

// add string s1 with string s2
// 文字列s1と文字列s2を追加する
char *addstr(const char *s1, const char *s2) {
    char *result = malloc(strlen(s1) + strlen(s2));
    // copy s1 to result
    strncpy(result, s1, strlen(s1));
    // copy s2 to result
    strncpy(result + strlen(s1), s2, strlen(s2));

    result[strlen(s1) + strlen(s2)] = '\0';

    return result;
}

// swap a with b　
void swap(int *a, int *b) {
    int tmp = *b;
    *b = *a;
    *a = tmp;
}

void ci_random(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]) {
    initialize(prec, m);
    int i, j, id_points[MAX_N], len_tour;

    for (i = 0; i < m; ++i) {
        tour[i] = prec[i];
    }
    tour[m] = tour[0];
    len_tour = m + 1;

    j = 0;
    for (i = 0; i < n; ++i) {
        if (flags[i] == 1) continue;
        id_points[j++] = i;
    }
    int size_id_points = n - m;

    shuffle(id_points, size_id_points, sizeof(int));

    for (i = 0; i < size_id_points; ++i) {
        double best_cost = INF;
        int best_pos = -1;
        for (j = 1; j < len_tour; ++j) {
            double cost_gain = dist(p[tour[j - 1]], p[id_points[i]]) + dist(p[id_points[i]], p[tour[j]]) - dist(p[tour[j - 1]], p[tour[j]]);
            if (cost_gain < best_cost) {
                best_cost = cost_gain;
                best_pos = j;
            }
        }
        insert(tour, &len_tour, best_pos, id_points[i]);
    }

    printf("Cheapest_insertion_random : %f\n", tour_length(p, n, tour));
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

void swap_points(struct point p[], int n, int tour[]) {
    tour[n] = tour[0];
    int index_points[n];
    int index_i, index_j;
    double delta;
    // show_array(tour, size_tour + 1);
    // printf("Distance tour start swap = %lf\n", tour_length(p, size_tour, tour));

    for (int i = 0; i < n; ++i) {
        index_points[i] = i;
    }

    shuffle(index_points, n, sizeof(int));

    for (int i = 0; i < n; ++i) {
        index_i = index_points[i];
        if (index_i == 0) continue;
        for (int j = i + 1; j < n; ++j) {
            index_j = index_points[j];
            if (index_j == 0 || index_j == index_i) continue;
            if (flags[tour[index_i]] || flags[tour[index_j]]) continue;

            if (index_i > index_j) swap(&index_i, &index_j);

            delta = 0.0;
            int pre_i = tour[index_i - 1];
            int cur_i = tour[index_i];
            int pos_i = tour[index_i + 1];
            int pre_j = tour[index_j - 1];
            int cur_j = tour[index_j];
            int pos_j = tour[index_j + 1];

            if (abs(index_i - index_j) == 1) {
                delta -= (dist(p[pre_i], p[cur_i]) + dist(p[cur_j], p[pos_j]));
                // printf("delta: %f", delta);
                delta += (dist(p[pre_i], p[cur_j]) + dist(p[cur_i], p[pos_j]));
                // printf("=> %f\n", delta);
            } else {
                delta -= (dist(p[pre_i], p[cur_i]) + dist(p[cur_i], p[pos_i]));
                delta -= (dist(p[pre_j], p[cur_j]) + dist(p[cur_j], p[pos_j]));
                // printf("delta: %f", delta);

                delta += (dist(p[pre_i], p[cur_j]) + dist(p[cur_j], p[pos_i]));
                delta += (dist(p[pre_j], p[cur_i]) + dist(p[cur_i], p[pos_j]));
                // printf("=> %f\n", delta);
            }
            // printf("swap: %d <-> %d : ", cur_i, cur_j);
            // printf("delta: %f\n", delta);

            if (delta < -EPSILON) {
                swap(&tour[index_i], &tour[index_j]);
                change = 1;
                double cur_distance = tour_length(p, n, tour);
                printf("Swap changed = %lf\n", cur_distance);
                if ((min_distance - cur_distance) > EPSILON) {
                    min_distance = cur_distance;
                    write_tour_data(tourFileName, n, tour);
                }
                return;
            }
        }
    }
}

void relocate(struct point p[], int n, int tour[]) {
    tour[n] = tour[0];
    int index_points[n];
    int index_i, index_j;
    double delta, cur_distance;
    // show_array(tour, size_tour + 1);
    // printf("Distance tour start swap = %lf\n", tour_length(p, size_tour, tour));

    for (int i = 0; i < n; ++i) {
        index_points[i] = i;
    }

    shuffle(index_points, n, sizeof(int));

    for (int i = 0; i < n; ++i) {
        index_i = index_points[i];
        if (index_i == 0) continue;
        if (flags[tour[index_i]]) continue;

        for (int j = 0; j < n; ++j) {
            index_j = index_points[j];
            if (index_j == index_i || index_j == (index_i - 1)) continue;

            delta = 0.0;
            int pre_i = tour[index_i - 1];
            int cur_i = tour[index_i];
            int pos_i = tour[index_i + 1];
            int cur_j = tour[index_j];
            int pos_j = tour[index_j + 1];

            delta = -dist(p[pre_i], p[cur_i]) - dist(p[cur_i], p[pos_i]) + dist(p[pre_i], p[pos_i]);
            delta += -dist(p[cur_j], p[pos_j]) + dist(p[cur_j], p[cur_i]) + dist(p[cur_i], p[pos_j]);

            if (delta < -EPSILON) {
                // printf("%d -> after %d - Delta: %lf\n", cur_i, cur_j, delta);
                int len = n + 1;
                delete (tour, &len, index_i);
                if (index_i < index_j) index_j--;
                insert(tour, &len, index_j + 1, cur_i);
                change = 1;
                cur_distance = tour_length(p, n, tour);
                printf("Relocate changed = %lf\n", cur_distance);
                // show_array(tour, size_tour + 1);
                return;
            }
        }
    }
    if ((min_distance - cur_distance) > EPSILON) {
        min_distance = cur_distance;
        write_tour_data(tourFileName, n, tour);
    }
}

void TwoOpt(struct point p[MAX_N], int n, int tour[MAX_N], int m, int prec[MAX_N]){
  int a,b,c,d;
  int i,j,k,l,g,h,r,z;
  int count;
  int countPrec=0,finish=0;
  double max_dist=0;
  int pc[MAX_N];
  double cur_distance;
  cur_distance = tour_length(p, n, tour);

  for(i=0;i<n;i++) //pcの初期化
    pc[i]=0;

  for(i=0;i<m;i++) //precの点を1としておく
    pc[prec[i]]=1;

  while(max_dist != (-1)*INF ){
    max_dist=(-1)*INF;

    for(i=0;i<=n-3;i++){
      j=i+1;
      count = 0;
      count += pc[tour[j]];

      for(k=i+2;k<=n-1;k++){
        l=(k+1)%n;

        count += pc[tour[k]];

        if(count > 1) break;

        a = tour[i];
        b = tour[j];
        c = tour[k];
        d = tour[l];

        if(dist(p[a],p[b])+dist(p[c],p[d])>dist(p[a],p[c])+dist(p[b],p[d])){
          if(max_dist < dist(p[a],p[b])+dist(p[c],p[d])-dist(p[a],p[c])-dist(p[b],p[d])){
            max_dist = dist(p[a],p[b])+dist(p[c],p[d])-dist(p[a],p[c])-dist(p[b],p[d]);
            g=j;
            h=k;
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

    double distance_tour = tour_length(p, n, tour);
    if ((cur_distance - distance_tour) > EPSILON) {
        printf("Changed by Two-Opt = %lf\n", distance_tour);
        cur_distance = distance_tour;
        change++;
        if ((min_distance - cur_distance) > EPSILON) {
            min_distance = cur_distance;
            write_tour_data(tourFileName, n, tour);
        }
        return;
    }
}

void TwoOpt_old(struct point p[MAX_N], int n, int tour[MAX_N], int m, int prec[MAX_N]){
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

void OrOpt1(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]){
  int a,b,c,d,e;
  int i0, i, i1, j, j1, g;
  int tmp, count=1;
  double k1, k2 ,k3, k4, k5, k6;
  double cur_distance;
  cur_distance = tour_length(p, n, tour);

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
      double distance_tour = tour_length(p, n, tour);
        if ((cur_distance - distance_tour) > EPSILON) {
            printf("Changed by orOpt1 = %lf\n", distance_tour);
            cur_distance = distance_tour;
            change++;
            if ((min_distance - cur_distance) > EPSILON) {
                min_distance = cur_distance;
                write_tour_data(tourFileName, n, tour);
            }
            return;
        }
  }
}

void OrOpt2(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]){
  int a,b,c,d,e,f;
  int i0, i, i1, i2, j, j1, g, tmp, tmp1;
  int count=1;
  double k1, k2 ,k3, k4, k5, k6, k7 , k8, k9, k10;
  double cur_distance;
  cur_distance = tour_length(p, n, tour);
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
        double distance_tour = tour_length(p, n, tour);
        if ((cur_distance - distance_tour) > EPSILON) {
            printf("Changed by orOpt2 = %lf\n", distance_tour);
            cur_distance = distance_tour;
            change++;
            if ((min_distance - cur_distance) > EPSILON) {
                min_distance = cur_distance;
                write_tour_data(tourFileName, n, tour);
            }
            return;
        }
  }
}


void OrOpt3(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]){
  int a,b,c,d,e,f,g;
  int i0, i, i1, i2, i3, j, j1, tmp, tmp1, tmp2, count=1;
  double k1, k2 ,k3, k4, k5, k6, k7 , k8, k9, k10, k11, k12;
  double min = INF;
  double cur_distance;
  cur_distance = tour_length(p, n, tour);

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
	    double distance_tour = tour_length(p, n, tour);
        if ((cur_distance - distance_tour) > EPSILON) {
            printf("Changed by orOpt3 = %lf\n", distance_tour);
            cur_distance = distance_tour;
            change++;
            if ((min_distance - cur_distance) > EPSILON) {
                min_distance = cur_distance;
                write_tour_data(tourFileName, n, tour);
            }
            return;
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
  if(S_tour > tour_length(p,n,tour)) {
	    change++;
		S_tour = tour_length(p,n,tour);
		sprintf(tourFileName,  "tour%08d.dat",++num);
		write_tour_data(tourFileName,n,tour);
		printf("Changed by orOpt %d = %lf\n", k, tour_length(p,n,tour));
  }
}


void local_search(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]) {
    change = 1;
    while (change) {
        if(tour_length(p,n,tour) < 100000){
            printf("Relocate!\n");
            relocate(p, n, tour);
            change = 0;
            printf("Swap!\n");
            swap_points(p, n, tour);
            printf("Or-Opt1!\n");
            OrOpt1(p, n, tour, m, prec);
            printf("Or-Opt2!\n");
            OrOpt2(p, n, tour, m, prec);
            printf("Or-Opt3!\n");
            OrOpt3(p, n, tour, m, prec);
            printf("2-Opt!\n");
            TwoOpt(p, n, tour, m, prec);
        } else if(tour_length(p,n,tour) > 500000){
            printf("Relocate!\n");
            relocate(p, n, tour);
            change = 0;
            printf("Swap!\n");
            swap_points(p, n, tour);
            printf("Or-Opt1!\n");
            OrOpt1(p, n, tour, m, prec);
            printf("Or-Opt2!\n");
            OrOpt2(p, n, tour, m, prec);
            printf("Or-Opt3!\n");
            OrOpt3(p, n, tour, m, prec);
            printf("2-Opt-Old!\n");
            TwoOpt_old(p, n, tour, m, prec);

        } else {
            printf("Relocate!\n");
            relocate(p, n, tour);
            change = 0;
            printf("Swap!\n");
            swap_points(p, n, tour);
            printf("Or-Opt3!\n");
            OrOpt3(p, n, tour, m, prec);
            printf("Or-Opt2!\n");
            OrOpt2(p, n, tour, m, prec);
            printf("Or-Opt1!\n");
            OrOpt1(p, n, tour, m, prec);
            printf("2-Opt-Old-!\n");
            TwoOpt_old(p, n, tour, m, prec);
        }

    }
}

void solve(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]) {
    double cur_distance, best_distance;

    ci_random(p,n,tour,m,prec);

    write_tour_data(tourFileName, n, tour);

    best_distance = tour_length(p, n, tour);
    min_distance = best_distance;

    for (int i = 0; i < max_i; ++i) {
        printf("i  %d :----------------\n", i);
        if (i >0)
            ci_random(p,n,tour,m,prec);
        local_search(p,n,tour,m,prec);
        

        cur_distance = tour_length(p, n, tour);
        if ((best_distance - cur_distance) > EPSILON) {
            printf("%d - Found improvement %lf\n", i, cur_distance);
            best_distance = cur_distance;
            min_distance = best_distance;
            write_tour_data(tourFileName, n, tour);
        }
    }

    printf("Min: %.2lf\n", min_distance);
    printf("Done Local Search!\n");
}

int main(int argc, char *argv[]) {
  int n;                   // 点の数 
  int m;                   // 順序制約に現れる点の数
  int k;
  struct point p[MAX_N];   // 各点の座標を表す配列 
  int tour[MAX_N];   // 巡回路を表現する配列
  int prec[MAX_N];   // 順序制約を表現する配列
  int i;

  if(argc != 2) {
    fprintf(stderr,"Usage: %s <tsp_filename>\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  // 点の数と各点の座標を1番目のコマンドライン引数で指定されたファイルから読み込む
  read_tsp_data(argv[1],p,&n,prec,&m);
/*
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
*/
   solve(p,n,tour,m,prec);
  
  exit(EXIT_SUCCESS);
}