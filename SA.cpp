
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>

#define T0 50000.0
#define T_end (1e-8)
#define q  0.98
#define L 1000


 double distance(double *city1,double *city2)
 {
     double x1 = *city1;
     double y1 = *(city1+1);
     double x2 = *(city2);
     double y2 = *(city2+1);
     double dis = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
     return dis;
 }


double path_len(int * arr, int N, double city_pos[][2])
{
     double path = 0;
     int index = *arr;
     for(int i=0;i<N-1;i++)
     {
         int index1 = *(arr+i);
         int index2 = *(arr+i+1);
         double dis = distance(city_pos[index1-1],city_pos[index2-1]);
         path += dis;
         //printf("%lf %lf %.2lf\n\n", city_pos[index1 - 1][0], city_pos[index1 - 1][1], dis);

     }
     int last_index = *(arr+N-1);
     int first_index = *arr;
     double last_dis = distance(city_pos[last_index-1],city_pos[first_index-1]);
     path = path + last_dis;
     return path;
 }

 void printList(int A[], int n){
    for(int i=0; i<n; i++){
        printf("%d ", A[i]);
    }
    printf("\n");
}

void init(int N, int city_list[])
{
     for(int i=0;i<N;i++){
         city_list[i] = i + 1;
     }
}

int isPossible(int city[], int constraint[], int nCity, int nConstraint){
    int cst[nCity];
    for(int i=0; i<nConstraint; i++) cst[constraint[i]] = 1;
    int cityCst[nConstraint];
    int index = 0;

    for(int i=0; i<nCity; i++){
       if(cst[city[i]] == 1){
        cityCst[index++] = city[i];
       }
    }

    for(int i=0; i<nConstraint; i++){
        if(cityCst[i] != constraint[i]){
            return 0;
        }
    }

    return 1;
}

void create_new(int N, int city_list[], int nConstraint, int constraint[])
{
  //  int cnt = 5;
   // do{
        double r1 = ((double)rand())/(RAND_MAX+1.0);
        double r2 = ((double)rand())/(RAND_MAX+1.0);

        int pos1 = (int)(N*r1);
        int pos2 = (int)(N*r2);
      //  printf("%d-%d\n", pos1, pos2);
        int temp = city_list[pos1];
        city_list[pos1] = city_list[pos2];
        city_list[pos2] = temp;

       // printf("%d\n", cnt);
    //}
    //while(!isPossible(city_list, constraint, N, nConstraint) && cnt > 0);

}




void initFirstResult(int result[], int N_CONSTRAINT, int constraint[], int N){

    int dd = (int) (N / N_CONSTRAINT) - 1;

    int dttr = rand() % dd + 1;
    int vis[N + 1];

    for(int i=0; i<=N; i++) vis[i] = 0;

    int tmp = 0;
    for(int i=1; i<N && tmp < N_CONSTRAINT; i+=dttr){
        result[i] = constraint[tmp];
        vis[constraint[tmp++]] = 1;
    }


    int trr[N - N_CONSTRAINT];

    int ind = 0;
    for(int i=0; i<N; i++){
        if(vis[i] != 1){
            trr[ind++] = i + 1;
        }
    }

    int curSize = N - N_CONSTRAINT;

    int index = 0, value = 1;


    while(index < N){
        if(result[index] > 0){
            index++;
            continue;
        }
        //printf("%d\n\n", curSize);
        int cIndex =  (rand() % curSize);


        result[index++] = trr[cIndex];

        for(int i=cIndex; i<curSize - 1; i++){
            trr[i] = trr[i + 1];
        }

        curSize--;

    }
}

int main(void)
{

   // INPUT
    FILE* f = fopen("inp.txt", "r");
    int N, N_CONSTRAINT;
    fscanf(f, "%d", &N);
    fscanf(f, "%d", &N_CONSTRAINT);
    int constraint[N_CONSTRAINT];
    for(int i=0; i<N_CONSTRAINT; i++){
        fscanf(f, "%d", &constraint[i]);
        constraint[i]++;
    }

    double city_pos[N][2];
    for(int i=0; i<N; i++){
        int x;
        fscanf(f, "%d", &x);
        fscanf(f, "%lf", &city_pos[i][0]);
        fscanf(f, "%lf", &city_pos[i][1]);
    }
    // END INPUT



    int city_list[N];


      srand((unsigned)time(NULL));
      time_t start,finish;
      start = clock();
      double T;
      int count = 0;
      T = T0;
      init(N, city_list);
      int city_list_copy[N];
      double f1,f2,df;
      double r;
      int result[N];

      for(int i=0; i<N; i++) result[i] = 0;


      initFirstResult(result, N_CONSTRAINT, constraint, N);
      memcpy(city_list,result,N*sizeof(int));


      while(T > T_end)
      {
          for(int i=0;i<L;i++)
          {
                memcpy(city_list_copy,city_list,N*sizeof(int));
                create_new(N, city_list, N_CONSTRAINT, constraint);
                f1 = path_len(city_list_copy, N, city_pos);
                f2 = path_len(city_list, N, city_pos);
                df = f2 - f1;

                if(df >= 0)
                {
                    r = ((double)rand())/(RAND_MAX);
                    if(exp(-df/T) <= r)
                    {

                        memcpy(city_list,city_list_copy,N*sizeof(int));
                        if(isPossible(city_list_copy, constraint, N, N_CONSTRAINT))
                            memcpy(result,city_list_copy,N*sizeof(int));
                    }
                }
         }
         T *= q;
         count++;
     }

     for(int i=0;i<N-1;i++)
     {
         printf("%d--->",result[i] - 1);
     }
     printf("%d\n",result[N-1]);
     double len = path_len(result, N, city_pos);
     printf("=> COST: %lf\n",len);

     fclose(f);
     return 0;
 }
