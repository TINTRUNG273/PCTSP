// GROUP 08: MEGANEs
// プログラムコンテストの第2ステージ

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// #include "utils.h"

#define MAX_N (int)1e4  // 点の数の最大値
#define INF (int)1e9    // 無限大の定義
#define EPSILON 1e-8    //ε 小さい正の値
#define SWAP(a, b)  \
    {               \
        int temp;   \
        temp = (a); \
        (a) = (b);  \
        (b) = temp; \
    }

char tourFileName[20];
int num = 0, max_iter = (int)1e9, s, change = INF, flags[MAX_N];
double min_distance = INF;

// random comparison function to create shuffle array
int rand_comparison(const void *a, const void *b) {
    (void)a;
    (void)b;
    return rand() % 2 ? +1 : -1;
}

// shuffle arrays with qsort
void shuffle(void *base, size_t nmemb, size_t size) {
    qsort(base, nmemb, size, rand_comparison);
}

// get substring from string s from position start to position end
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

typedef struct point {
    int x;
    int y;
} point;

double dist(point p, point q) {  // pとq の間の距離を計算
    return sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y));
}

double tour_length(point p[], int n, int tour[]) {
    int i;
    double sum = 0.0;
    for (i = 0; i < n; i++) sum += dist(p[tour[i]], p[tour[(i + 1) % n]]);
    return sum;  // 総距離が関数の戻り値
}

// initialize the values of variables
void initialize(int prec[], int number_prec) {
    for (int i = 0; i < number_prec; ++i) {
        int id_node = prec[i];
        flags[id_node] = 1;
    }
    sprintf(tourFileName, "tour%08d.dat", 1);
}

// insert value at position k in array
void insert(int tour[], int *len, int k, int value) {
    int i;

    if (k < 0 || k > *len) {
        printf("Error in insert: out of range\n");
    }

    for (i = *len; i > k; i--) {
        tour[i] = tour[i - 1];
    }
    tour[k] = value;
    (*len)++;
}

// delete element at position k in array
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

void reverse(int arr[], int start, int end) {
    for (int i = 0; i < (end - start + 1) / 2; ++i) {
        swap(&arr[start + i], &arr[end - i]);
    }
}

void show_array(int a[], int len) {
    int i;
    printf("[ ");
    for (i = 0; i < len - 1; i++) {
        printf("%d ", a[i]);
    }
    printf("%d ]\n", a[i]);
}

void read_tsp_data(char *filename, point p[], int *np, int prec[], int *mp) {
    FILE *fp;
    char buff[500];
    int i;

    if ((fp = fopen(filename, "rt")) == NULL) {  // 指定ファイルを読み込み用に開く
        fprintf(stderr, "Error: File %s open failed.\n", filename);
        exit(EXIT_FAILURE);
    }

    while ((fgets(buff, sizeof(buff), fp) != NULL)  // PRECEDENCE_CONSTRAINTS:で始まる行に出会う
           && (strncmp("PRECEDENCE_CONSTRAINTS:", buff, 23) != 0))
        ;  // まで読み飛ばす.
    if (strncmp("PRECEDENCE_CONSTRAINTS:", buff, 23) == 0) {
        sscanf(buff + 24, "%d", mp);
        for (i = 0; i < *mp; i++) {
            fscanf(fp, "%d ", &prec[i]);
        }

    } else {
        fprintf(stderr, "Error: There is no precedence constraint in file %s.\n", filename);
        exit(EXIT_FAILURE);
    }

    while ((fgets(buff, sizeof(buff), fp) != NULL)  // DIMENSION で始まる行に出会う
           && (strncmp("DIMENSION", buff, 9) != 0))
        ;                               // まで読み飛ばす.
    sscanf(buff, "DIMENSION: %d", np);  // 点の数 *np を読み込む

    while ((fgets(buff, sizeof(buff), fp) != NULL)  // NODE_COORD_SECTIONで始まる
           && (strncmp("NODE_COORD_SECTION", buff, 18) != 0))
        ;                        // 行に出会うまで,
                                 // 読み飛ばす.
    for (i = 0; i < *np; i++) {  // i=0 から i=(*np)-1まで
        if (fgets(buff, sizeof(buff), fp) != NULL)
            sscanf(buff, "%*d %d %d", &(p[i].x), &(p[i].y));  // i番目の点の座標を読み込む
    }

    fclose(fp);
}

void write_tour_data(char *filename, int n, int tour[]) {
    FILE *fp;
    int i;

    // 構築した巡回路をfilenameという名前のファイルに書き出すためにopen
    if ((fp = fopen(filename, "wt")) == NULL) {
        fprintf(stderr, "Error: File %s open failed.\n", filename);
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "%d\n", n);
    for (i = 0; i < n; i++) {
        fprintf(fp, "%d ", tour[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
}

int check_prec(int prec[], int r, int m) {
    // int i;
    // for (i = 0; i < m; i++)
    //     if (prec[i] == r) {
    //         return 1;
    //     }
    // return 0;
    return flags[r];
}

// choose the position with the shortest distance to insert
void cheapest_heuristic_random(point p[], int prec[], int tour[], int number_prec, int size_tour) {
    initialize(prec, number_prec);
    int i, j, id_points[MAX_N], len_tour;

    // create an initial tour containing PRECEDENCE CONSTRAINTS points
    for (i = 0; i < number_prec; ++i) {
        tour[i] = prec[i];
    }
    tour[number_prec] = tour[0];
    len_tour = number_prec + 1;

    // save points that are not part of PRECEDENCE CONSTRAINTS to id_points
    j = 0;
    for (i = 0; i < size_tour; ++i) {
        if (flags[i] == 1) continue;
        id_points[j++] = i;
    }
    int size_id_points = size_tour - number_prec;

    // shuffling array id_points
    shuffle(id_points, size_id_points, sizeof(int));

    // insert each point in the id_points array into the tour at the position with the least change in distance
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

    // check(tour, size_tour, prec, number_prec);
    printf("Distance tour cheapest_heuristic_random : %f\n", tour_length(p, size_tour, tour));
}

// swap each pair of points in the tour
void swap_points(int tour[], point p[], int size_tour) {
    tour[size_tour] = tour[0];
    int index_points[size_tour];
    int index_i, index_j;
    double delta, cur_distance;
    // show_array(tour, size_tour + 1);
    // printf("Distance tour start swap = %lf\n", tour_length(p, size_tour, tour));
    cur_distance = tour_length(p, size_tour, tour);

    for (int i = 0; i < size_tour; ++i) {
        index_points[i] = i;
    }

here:
    shuffle(index_points, size_tour, sizeof(int));

    for (int i = 0; i < size_tour; ++i) {
        index_i = index_points[i];
        if (index_i == 0) continue;
        for (int j = i + 1; j < size_tour; ++j) {
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

            if (abs(index_i - index_j) == 1) {  // case of 2 adjacent points
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

            if (delta < -EPSILON) {  // if there is an improvement in the distance for the tour
                swap(&tour[index_i], &tour[index_j]);
                change = 1;
                // cur_distance = tour_length(p, size_tour, tour);
                // printf("Swap changed = %lf\n", cur_distance);
                // if ((min_distance - cur_distance) > EPSILON) {
                //     min_distance = cur_distance;
                //     write_tour_data(tourFileName, size_tour, tour);
                // }
                goto here;
                return;
            }
        }
    }

    double distance = tour_length(p, size_tour, tour);
    if (distance + EPSILON < cur_distance) {
        printf("Swap changed = %lf\n", distance);
    }
}

// change the position of each point in the tour
void relocate(int tour[], point p[], int size_tour) {
    tour[size_tour] = tour[0];
    int index_points[size_tour];
    int index_i, index_j;
    double delta, cur_distance;
    // show_array(tour, size_tour + 1);
    // printf("Distance tour start swap = %lf\n", tour_length(p, size_tour, tour));
    cur_distance = tour_length(p, size_tour, tour);

    for (int i = 0; i < size_tour; ++i) {
        index_points[i] = i;
    }

here:
    shuffle(index_points, size_tour, sizeof(int));
    // printf("Shuffle\n");

    for (int i = 0; i < size_tour; ++i) {
        index_i = index_points[i];
        // printf("index i: %d\n", index_i);
        if (index_i == 0) continue;
        if (flags[tour[index_i]]) continue;

        for (int j = 0; j < size_tour; ++j) {
            index_j = index_points[j];
            if (index_j == index_i || index_j == (index_i - 1)) continue;  // The new location and the old location are the same
            // printf("index j: %d\n", index_j);

            delta = 0.0;
            int pre_i = tour[index_i - 1];
            int cur_i = tour[index_i];
            int pos_i = tour[index_i + 1];
            int cur_j = tour[index_j];
            int pos_j = tour[index_j + 1];

            delta = -dist(p[pre_i], p[cur_i]) - dist(p[cur_i], p[pos_i]) + dist(p[pre_i], p[pos_i]);
            delta += -dist(p[cur_j], p[pos_j]) + dist(p[cur_j], p[cur_i]) + dist(p[cur_i], p[pos_j]);

            // printf("delta:%lf\n", delta);

            if (delta < -EPSILON) {  // if there is an improvement in the distance for the tour
                // printf("%d -> after %d - Delta: %lf\n", cur_i, cur_j, delta);
                int len = size_tour + 1;
                delete (tour, &len, index_i);
                if (index_i < index_j) index_j--;
                insert(tour, &len, index_j + 1, cur_i);
                change = 1;
                // cur_distance = tour_length(p, size_tour, tour);
                // printf("Relocate changed = %lf\n", cur_distance);
                // show_array(tour, size_tour + 1);

                goto here;
                return;
            }
        }
    }
    double distance = tour_length(p, size_tour, tour);
    if (distance + EPSILON < cur_distance) {
        printf("Relocate changed = %lf\n", distance);
    }
    // if ((min_distance - cur_distance) > EPSILON) {
    //     min_distance = cur_distance;
    //     write_tour_data(tourFileName, size_tour, tour);
    // }
}

void TwoOpt(point p[], int size_tour, int tour[]) {
    tour[size_tour] = tour[0];
    int index_points[size_tour];
    int index_i, index_j, count_pre;
    double delta, cur_distance;
    // show_array(tour, size_tour + 1);
    // printf("Distance tour start Two-opt = %lf\n", tour_length(p, size_tour, tour));
    cur_distance = tour_length(p, size_tour, tour);

    for (int i = 0; i < size_tour; ++i) {
        index_points[i] = i;
    }

here:
    count_pre = 0;
    shuffle(index_points, size_tour, sizeof(int));

    for (int i = 0; i < size_tour; ++i) {
        index_i = index_points[i];
        if (index_i == 0) continue;
        int pre_i = tour[index_i - 1];
        int cur_i = tour[index_i];
        count_pre = flags[cur_i];

        for (index_j = index_i + 1; index_j < size_tour; ++index_j) {
            int cur_j = tour[index_j];
            int pos_j = tour[index_j + 1];
            count_pre += flags[cur_j];
            if (count_pre >= 2) break;

            delta = 0;
            delta -= (dist(p[pre_i], p[cur_i]) + dist(p[cur_j], p[pos_j]));
            delta += (dist(p[pre_i], p[cur_j]) + dist(p[cur_i], p[pos_j]));

            if (delta < -EPSILON) {  // if there is an improvement in the distance for the tour
                reverse(tour, index_i, index_j);
                change = 1;
                // cur_distance = tour_length(p, size_tour, tour);

                // printf("2-opt changed = %lf\n", cur_distance);
                // if ((min_distance - cur_distance) > EPSILON) {
                //     min_distance = cur_distance;
                //     write_tour_data(tourFileName, size_tour, tour);
                // }
                goto here;
                return;
            }
        }
    }
    // if ((min_distance - cur_distance) > EPSILON) {
    //     min_distance = cur_distance;
    //     write_tour_data(tourFileName, size_tour, tour);
    // }
    double distance = tour_length(p, size_tour, tour);
    if (distance + EPSILON < cur_distance) {
        printf("2-opt changed = %lf\n", distance);
    }
}

void OrOpt1(struct point p[MAX_N], int n, int tour[MAX_N], int m, int prec[MAX_N]) {
    int a, b, c, d, e;
    int i0, i, i1, j, j1, g, tmp;
    double k1, k2, k3, k4, k5, k6;
    double cur_distance;

here:
    cur_distance = tour_length(p, n, tour);
    int count = 1;

    while (count > 0) {
        count = 0;
        for (i = 0; i < n; i++) {
            count = 0;
            i0 = i - 1;
            if (i0 < 0) i0 = n - 1;
            i1 = (i + 1) % n;
            if (check_prec(prec, tour[i], m) > 0) continue;
            for (j = 0; j < n; j++) {
                if (check_prec(prec, tour[i], m) > 0) continue;
                j1 = (j + 1) % n;
                a = tour[i0];
                b = tour[i1];
                c = tour[i];
                d = tour[j];
                e = tour[j1];
                if (j != i && j1 != i) {
                    k1 = dist(p[a], p[c]);
                    k2 = dist(p[c], p[b]);
                    k3 = dist(p[d], p[e]);
                    k4 = dist(p[a], p[b]);
                    k5 = dist(p[d], p[c]);
                    k6 = dist(p[c], p[e]);
                    if (k1 + k2 + k3 > k4 + k5 + k6) {
                        count++;
                        tmp = tour[i];
                        if (i < j) {
                            for (g = i; g < j; g++) tour[g] = tour[g + 1];
                            tour[j] = tmp;
                            i--;
                        } else {
                            for (g = i; g > j1; g--) tour[g % n] = tour[(g - 1) % n];
                            tour[j1] = tmp;
                        }
                        break;
                    }
                }
            }
        }

        double distance_tour = tour_length(p, n, tour);
        if ((cur_distance - distance_tour) > EPSILON) {
            // check(tour, n, prec, m);

            printf("Changed by orOpt1 = %lf\n", distance_tour);
            cur_distance = distance_tour;
            change++;
            // if ((min_distance - cur_distance) > EPSILON) {
            //     min_distance = cur_distance;
            //     write_tour_data(tourFileName, n, tour);
            // }
            goto here;
            return;
        }
    }
}

void OrOpt2(struct point p[MAX_N], int n, int tour[MAX_N], int m, int prec[MAX_N]) {
    int a, b, c, d, e, f;
    int i0, i, i1, i2, j, j1, g, tmp, tmp1;
    double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10;
    double cur_distance;

here:
    cur_distance = tour_length(p, n, tour);
    int count = 1;

    while (count > 0) {
        count = 0;
        for (i = 0; i < n; i++) {
            if (check_prec(prec, tour[i], m) > 0 || check_prec(prec, tour[i + 1], m) > 0) {
                continue;
            }
            i0 = i - 1;
            if (i0 < 0) i0 = n - 1;
            i1 = (i + 1) % n;
            i2 = (i1 + 1) % n;

            for (j = i + 2; j < n; j++) {
                if (check_prec(prec, tour[i], m) > 0 || check_prec(prec, tour[i + 1], m) > 0) continue;
                j1 = (j + 1) % n;
                a = tour[i0];
                b = tour[i2];
                c = tour[i];
                d = tour[i1];
                e = tour[j];
                f = tour[j1];
                if (j != i && j1 != i) {
                    k1 = dist(p[a], p[c]);
                    k2 = dist(p[d], p[b]);
                    k3 = dist(p[e], p[f]);
                    k4 = dist(p[a], p[b]);
                    k5 = dist(p[e], p[c]);
                    k6 = dist(p[d], p[f]);
                    k7 = dist(p[e], p[d]);
                    k8 = dist(p[c], p[f]);
                    k9 = dist(p[a], p[d]);
                    k10 = dist(p[c], p[b]);

                    if (k1 + k2 + k3 > k4 + k5 + k6) {
                        count++;
                        tmp = tour[i];
                        tmp1 = tour[(i + 1) % n];
                        for (g = i; g < j - 1; g++) tour[g] = tour[(g + 2) % n];
                        if (k5 + k6 > k7 + k8) {
                            tour[j - 1] = tmp1;
                            tour[j] = tmp;
                        } else {
                            tour[j - 1] = tmp;
                            tour[j] = tmp1;
                        }
                    }
                }
            }
        }

        double distance_tour = tour_length(p, n, tour);
        if ((cur_distance - distance_tour) > EPSILON) {
            // check(tour, n, prec, m);

            printf("Changed by orOpt2 = %lf\n", distance_tour);
            cur_distance = distance_tour;
            change++;
            // if ((min_distance - cur_distance) > EPSILON) {
            //     min_distance = cur_distance;
            //     write_tour_data(tourFileName, n, tour);
            // }
            goto here;
            return;
        }
    }
}

void OrOpt3(struct point p[MAX_N], int n, int tour[MAX_N], int m, int prec[MAX_N]) {
    int a, b, c, d, e, f, g;
    int i0, i, i1, i2, i3, j, j1, tmp, tmp1, tmp2;
    double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12;
    double cur_distance;

here:
    cur_distance = tour_length(p, n, tour);
    int count = 1;
    double min = INF;

    while (count > 0) {
        count = 0;
        for (i = 0; i < n; i++) {
            if (check_prec(prec, tour[i], m) > 0 || check_prec(prec, tour[i + 1], m) > 0) {
                continue;
            }
            i0 = i - 1;
            if (i0 < 0) i0 = n - 1;
            i1 = (i + 1) % n;
            i2 = (i1 + 1) % n;
            i3 = (i1 + 2) % n;

            for (j = i + 3; j < n; j++) {
                if ((check_prec(prec, tour[i], m) > 0 || check_prec(prec, tour[i + 1], m) > 0) || check_prec(prec, tour[i + 2], m) > 0) continue;
                j1 = (j + 1) % n;
                a = tour[i0];
                b = tour[i3];
                c = tour[i];
                d = tour[i1];
                e = tour[i2];
                f = tour[j];
                g = tour[j1];
                if (j != i && j1 != i) {
                    k1 = dist(p[a], p[c]);
                    k2 = dist(p[e], p[b]);
                    k3 = dist(p[f], p[g]);
                    k4 = dist(p[a], p[b]);
                    k5 = dist(p[f], p[c]);
                    k6 = dist(p[e], p[g]);
                    k7 = dist(p[f], p[c]);
                    k8 = dist(p[f], p[d]);
                    k9 = dist(p[f], p[e]);
                    k10 = dist(p[g], p[e]);
                    k11 = dist(p[g], p[d]);
                    k12 = dist(p[g], p[e]);

                    if (k1 + k2 + k3 > k4 + k5 + k6) {
                        count++;
                        tmp = tour[i];
                        tmp1 = tour[(i + 1) % n];
                        tmp2 = tour[(i + 2) % n];
                        for (g = i; g < j - 2; g++) tour[g] = tour[(g + 3) % n];
                        min = INF;

                        if (min > k7 + k12) {
                            min = k7 + k12;
                            tour[j - 2] = tmp;
                            tour[j - 1] = tmp1;
                            tour[j] = tmp2;
                        }

                        if (min > k8 + k12) {
                            min = k8 + k12;
                            tour[j - 2] = tmp1;
                            tour[j - 1] = tmp;
                            tour[j] = tmp2;
                        }
                        if (min > k7 + k11) {
                            min = k7 + k11;
                            tour[j - 2] = tmp;
                            tour[j - 1] = tmp2;
                            tour[j] = tmp1;
                        }
                        if (min > k8 + k10) {
                            min = k8 + k10;
                            tour[j - 2] = tmp1;
                            tour[j - 1] = tmp2;
                            tour[j] = tmp;
                        }
                        if (min > k9 + k11) {
                            min = k9 + k11;
                            tour[j - 2] = tmp2;
                            tour[j - 1] = tmp;
                            tour[j] = tmp1;
                        }
                        if (min > k9 + k10) {
                            min = k9 + k10;
                            tour[j - 2] = tmp2;
                            tour[j - 1] = tmp1;
                            tour[j] = tmp;
                        }
                    }
                }
            }
        }

        double distance_tour = tour_length(p, n, tour);
        if ((cur_distance - distance_tour) > EPSILON) {
            // check(tour, n, prec, m);

            printf("Changed by orOpt3 = %lf\n", distance_tour);
            cur_distance = distance_tour;
            change++;
            // if ((min_distance - cur_distance) > EPSILON) {
            //     min_distance = cur_distance;
            //     write_tour_data(tourFileName, n, tour);
            // }
            goto here;
            return;
        }
    }
}

void call(int f_i, point p[], int tour[], int size_tour, int prec[], int number_prec) {
    switch (f_i) {
        case 0:
            printf("2-Opt!\n");
            TwoOpt(p, size_tour, tour);
            break;
        case 1:
            printf("Swap!\n");
            swap_points(tour, p, size_tour);
            break;
        case 2:
            printf("Or-Opt1!\n");
            OrOpt1(p, size_tour, tour, number_prec, prec);
            break;
        case 3:
            printf("Or-Opt2!\n");
            OrOpt2(p, size_tour, tour, number_prec, prec);
            break;
        case 4:
            printf("Or-Opt3!\n");
            OrOpt3(p, size_tour, tour, number_prec, prec);
            break;
        default:
            printf("Relocate!\n");
            relocate(tour, p, size_tour);
            break;
    }
}

// local minimum value search
void local_search(point p[], int tour[], int size_tour, int prec[], int number_prec) {
    change = 1;
    int choose[6] = {0, 1, 2, 3, 4, 5};
    while (change) {
        shuffle(choose, 6, sizeof(int));
        change = 0;
        for (int i = 0; i < 6; i++) {
            call(choose[i], p, tour, size_tour, prec, number_prec);
        }
    }
}

// solver function
void solve(point p[], int tour[], int size_tour, int prec[], int number_prec) {
    double cur_distance, best_distance;

    cheapest_heuristic_random(p, prec, tour, number_prec, size_tour);

    write_tour_data(tourFileName, size_tour, tour);

    best_distance = tour_length(p, size_tour, tour);
    // min_distance = best_distance;

    for (int iter = 0; iter < max_iter; ++iter) {
        printf("iter %d:----------------\n", iter);
        if (iter > 0)
            cheapest_heuristic_random(p, prec, tour, number_prec, size_tour);

        local_search(p, tour, size_tour, prec, number_prec);
        // check(tour, size_tour, prec, number_prec);

        cur_distance = tour_length(p, size_tour, tour);
        if ((best_distance - cur_distance) > EPSILON) {
            printf("%d - Found improvement %.2lf\n", iter, cur_distance);
            best_distance = cur_distance;
            // min_distance = best_distance;
            write_tour_data(tourFileName, size_tour, tour);
        }
    }

    printf("Min: %.2lf\n", best_distance);
    printf("Done Local Search!\n");
}

int main(int argc, char *argv[]) {
    int n;            // 点の数
    int m;            // 順序制約に現れる点の数
    point p[MAX_N];   // 各点の座標を表す配列
    int tour[MAX_N];  // 巡回路を表現する配列
    int prec[MAX_N];  // 順序制約を表現する配列
    int i, j;
    srand(0);

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <tsp_filename> (max_iter)\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // ----------------------------------------------------------------
    // if (argc >= 3) {
    //     max_iter = 100;
    //     max_iter = atoi(argv[2]);
    // }
    //-----------------------------------------------------------------

    // 点の数と各点の座標を1番目のコマンドライン引数で指定されたファイルから読み込む
    read_tsp_data(argv[1], p, &n, prec, &m);

    solve(p, tour, n, prec, m);

    exit(EXIT_SUCCESS);
    return 0;
}
