#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>
#define T0 50000.0  // 初始温度
#define T_end (1e-8)
#define q  0.98   // 退火系数
#define L 1000  // 每个温度时的迭代次数，即链长
//#define N 31  // 城市数量
//int city_list[N]; // 用于存放一个解

// 中国31个城市坐标

 //函数声明
double distance(double *,double *); // 计算两个城市距离
double path_len(int *, double **, int);  // 计算路径长度
void  init(int *, int);  //初始化函数
void create_new(int *, int); // 产生新解
// 距离函数
double distance(double * city1,double * city2)
{
    double x1 = *city1;
    double y1 = *(city1+1);
    double x2 = *(city2);
    double y2 = *(city2+1);
    double dis = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    return dis;
}

// 计算路径长度
double path_len(int * arr, double city_pos[][2], int N)
{
    double path = 0; // 初始化路径长度
    int index = *arr; // 定位到第一个数字(城市序号)
    for(int i=0;i<N-1;i++)
    {
        int index1 = *(arr+i);
        int index2 = *(arr+i+1);
        double dis = distance(city_pos[index1-1],
                                city_pos[index2-1]);

         path += dis;
     }
    int last_index = *(arr+N-1); // 最后一个城市序号
    int first_index = *arr; // 第一个城市序号
    double last_dis = distance(city_pos[last_index-1],
                                city_pos[first_index-1]);
    path = path + last_dis;
    return path; // 返回总的路径长度
}

// 初始化函数
void init(int city_list[], int N)
{
    for(int i=0;i<N;i++)
        city_list[i] = i+1;  // 初始化一个解
}

// 产生一个新解
// 此处采用随机交换两个位置的方式产生新的解
void create_new(int city_list[], int N)
{
    double r1 = ((double)rand())/(RAND_MAX+1.0);
    double r2 = ((double)rand())/(RAND_MAX+1.0);
    int pos1 = (int)(N*r1); //第一个交叉点的位置
    int pos2 = (int)(N*r2);
    int temp = city_list[pos1];
    city_list[pos1] = city_list[pos2];
    city_list[pos2] = temp;   // 交换两个点
}



// 主函数
int main(void)
{
    srand((unsigned)time(NULL)); //初始化随机数种子
    time_t start,finish;
    start = clock(); // 程序运行开始计时
    double T;
    int count = 0; // 记录降温次数
    T = T0; //初始温度
    //read file
    int N;
    FILE* fp;
    fp = fopen("../TSP_data/wi29.tsp","rb");
    double tmp;
    fscanf(fp,"%d",&N);
    double city_pos[N][2];
    for(int i=0; i<N; i++)
  	{
  		for(int j=0; j<3; j++)
  		{
        if (j==0) fscanf(fp,"%lf",&tmp);
        else fscanf(fp,"%lf", &city_pos[i][j-1]);
  		}
  	}
    //
    int city_list[N], city_list_copy[N]; // 用于保存原始解
    init(city_list, N); //初始化一个解
    double f1,f2,df; //f1为初始解目标函数值，
                     //f2为新解目标函数值，df为二者差值


    double r; // 0-1之间的随机数，用来决定是否接受新解
    while(T > T_end) // 当温度低于结束温度时，退火结束
    {
        for(int i=0;i<L;i++)
        {
            // 复制数组
            memcpy(city_list_copy,city_list,N*sizeof(int));
            create_new(city_list, N); // 产生新解
            f1 = path_len(city_list_copy, city_pos, N);
            f2 = path_len(city_list, city_pos, N);
            df = f2 - f1;
            // 以下是Metropolis准则
            if(df >= 0)
            {
                r = ((double)rand())/(RAND_MAX);
                if(exp(-df/T) <= r) // 保留原来的解
                {
                    memcpy(city_list,city_list_copy,N*sizeof(int));
                }
            }
        }
        T *= q; // 降温
        count++;
    }
    finish = clock(); // 退火过程结束
    double duration = ((double)(finish-start))/CLOCKS_PER_SEC; // 计算时间
    printf("模拟退火算法，初始温度T0=%.2f,降温系数q=%.2f,每个温度迭代%d次,共降温%d次，得到的TSP最优路径为:\n",T0,q,L,count);
    for(int i=0;i<N-1;i++)  // 输出最优路径
    {
        printf("%d--->",city_list[i]);
    }
    printf("%d\n",city_list[N-1]);
    double len = path_len(city_list, city_pos, N); // 最优路径长度
    printf("最优路径长度为:%lf\n",len);
    printf("程序运行耗时:%lf秒.\n",duration);
    return 0;
}
