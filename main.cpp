#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <time.h>
#include <fstream>
#include <stack>
#include <math.h>

using namespace std;

typedef pair<int,int> POINT;

vector<POINT> point;

bool cmp_descend(POINT a,POINT b);
bool cmp_y_descend(POINT &a, POINT &b);
double multi(POINT &p1,POINT &p2,POINT &p0);
double Distance(POINT p1,POINT p2);

void data_generalization(vector<pair<int,int> > &point,int n,int range){
    for(int i = 0; i < n; i++)
    {
        int x = rand() % range;
        int y = rand() % range;
        point.push_back(make_pair(x,y));
        for(int j = 0; j < i;j++){
            if((x == point[j].first) && (y == point[j].second)){
                point.pop_back();i--;
                break;
            }
        }

    }
}

void CH(vector<POINT> &Convex_Hull_point,string file){
    ofstream out2(file.c_str());
    int n=Convex_Hull_point.size();
    sort(Convex_Hull_point.begin(),Convex_Hull_point.end());
    POINT start=Convex_Hull_point[0],Mid=Convex_Hull_point[n-1];
    vector<POINT> S_U,S_L;
    for (int i=1;i<n-1;i++){
        if (multi(Convex_Hull_point[i],start,Mid)>0) S_U.push_back(Convex_Hull_point[i]);
        else S_L.push_back(Convex_Hull_point[i]);
    }
    sort(S_U.begin(),S_U.end());
    sort(S_L.begin(),S_L.end(),cmp_descend);
    out2<<start.first<<" "<<start.second<<endl;
    for (unsigned int i=0;i<S_U.size();i++)
        out2<<S_U[i].first<<" "<<S_U[i].second<<endl;
    out2<<Mid.first<<" "<<Mid.second<<endl;
    for (unsigned int i=0;i<S_L.size();i++)
        out2<<S_L[i].first<<" "<<S_L[i].second<<endl;
}

vector<POINT> BruteForce(vector<POINT> &point)
{
    int n=point.size();
    int flag[n];
    //将所有点标志位设1,表示没有点被排除
    fill(flag,flag+n,1);
    //找纵坐标最小的点,若纵坐标相等找横坐标最小的点
    int u = 0;
    for(int i = 1; i < n; i++)
        if((point[i].second<point[u].second)||(point[i].second==point[u].second&&point[i].first<point[u].first))
            u = i;

    for(int j = 0; j < n; j++){
        if(j != u){
            if(flag[j] == 0)
                continue;
            for(int k=0;k<n;k++){
                if((k!=u)&&(k!=j)){
                    if(flag[k]==0)
                        continue;
                    for(int l=0;l<n;l++){
                        if((l!=u)&& (l!=j)&&(l!=k)){
                            if(flag[l] == 0)
                                continue;
                            //对于任意一条边，判断另两个点是否在同一侧，如果都在同一侧则删除该节点
                            vector<vector<int> > operator_flag(3,vector<int>(3,0));
                            operator_flag[0][1]=multi(point[l],point[j],point[u]);
                            operator_flag[0][2]=multi(point[k],point[j],point[u]);
                            if(operator_flag[0][1]*operator_flag[0][2] >= 0)
                                operator_flag[0][0]=1;

                            operator_flag[1][1]=multi(point[l],point[k],point[j]);
                            operator_flag[1][2]=multi(point[u],point[k],point[j]);
                            if(operator_flag[1][1]*operator_flag[1][2] >= 0)
                                operator_flag[1][0]=1;

                            operator_flag[2][1]=multi(point[l],point[u],point[k]);
                            operator_flag[2][2]=multi(point[j],point[u],point[k]);
                            if(operator_flag[2][1]*operator_flag[2][2] >= 0)
                                operator_flag[2][0]=1;

                            if(operator_flag[0][0]&operator_flag[1][0]&operator_flag[2][0])
                                flag[l] = 0;
                        }
                    }
                }
            }
        }
    }
    vector<POINT> convex;
    for (int i=0;i<n;i++){
        if (flag[i]){
            convex.push_back(point[i]);
        }
    }
    return convex;
}

void ascend_sort(vector<POINT> &point,vector<double> &theta){
    int n=point.size();
    for (int i=1;i<n;i++){
        double temp=atan2(point[i].second,point[i].first);
        theta[i]=temp<0?temp+M_PI:temp;
    }
    for (int i=0;i<n-1;i++){//极角排序,一定要把原点加入
        for (int j=i+1;j<n;j++){
            if (theta[i]>theta[j]){
               swap(point[i],point[j]);
               swap(theta[i],theta[j]);
            }
            else{
                if (theta[i]==theta[j]){
                    if (Distance(point[i],point[0])<=Distance(point[j],point[0])) swap(point[i],point[j]);
                }
            }
        }
    }
}

void sort_point_by_angle(vector<POINT> &point){
    int u = 0,n=point.size();
    for(int i=0; i < n; i++)//找出纵坐标最小的点，如果纵坐标相等找横坐标最小的点
        if((point[i].second < point[u].second) || (point[i].second == point[u].second && point[i].first < point[u].first))
            u = i;
    swap(point[0],point[u]);
    for(int i = 1; i < n; i++) {//按照极角从小到大将点进行排序
        point[i].first  = point[i].first  - point[0].first;
        point[i].second = point[i].second - point[0].second;
    }
    vector<double> theta(n,-1);
    ascend_sort(point,theta);
    //ofstream out1("sorted_point.txt");
    for(int i = 1; i < n; i++){
        point[i].first  = point[i].first  + point[0].first;
        point[i].second = point[i].second + point[0].second;
        //out1<<"("<<point[i].first <<","<<point[i].second<<")"<<theta[i]<<endl;
    }
}

vector<POINT> GrahamScan(vector<POINT> &point)
{
    int i;
    int top = 0;
    int n=point.size();
    vector<POINT> convex;
    int Stack[n];
    fill(Stack,Stack+n,0);
    sort_point_by_angle(point);
    //cout<<point[0].first<<" "<<point[0].second<<endl;
    for(i = 0; i <= 2;i++)//用栈来存储初始的三个点
        Stack[i] = i;
    top = 2;

    for(i = 3; i < n; i++)//不断判断栈顶两个元素与当前元素是否构成非左旋转，如果构成非左，则删除栈顶元素继续判断
    {
        while(multi(point[i], point[Stack[top]], point[Stack[top-1]]) > 0)
        {
            if(top == 0)
                break;
            top--;
        }
        if (multi(point[i], point[Stack[top]], point[Stack[top-1]]) == 0) top--;
        top++;
        Stack[top] = i;
    }

    while(multi(point[0], point[Stack[top]], point[Stack[top-1]]) > 0)
    {
        if(top == 0)
            break;
        top--;
    }
    if (multi(point[0], point[Stack[top]], point[Stack[top-1]]) == 0) top--;

    for(int i = 0; i <= top; i++){
        convex.push_back(point[Stack[i]]);
    }
    return convex;
}

vector<POINT> DivideConvex(vector<POINT> &point)
{
    int i;
    int n=point.size();

    int mid_n = n/2;
    vector<POINT> temp1, temp2, convex;

    if(n <= 3){
        return point;
    } //如果点集就剩三个元素，返回这三个元素
    sort(point.begin(),point.end(),cmp_y_descend);

    vector<POINT> left;
    vector<POINT> right;

    for(i = 0; i < mid_n; i++) //将点按照横坐标分成均等的左右两个集合
        left.push_back(point[i]);
    for(i = mid_n; i < n; i++)
        right.push_back(point[i]);

    temp1 = DivideConvex(left); //递归地在左右两个集合上进行处理
    temp2 = DivideConvex(right);

    for(unsigned int i = 0; i < temp2.size(); i++)
        temp1.push_back(temp2[i]); //将两个集合合并
    convex = GrahamScan(temp1); //在合并后的集合上应用GrahamScan算法

    return convex;
}

void Solve_solution_CH(vector<POINT> (*func)(vector<POINT> &point),vector<POINT> &point,string file){
    int n=point.size();
    //POINT right_top=point[n-1];
    //cout<<right_top.first<<" "<<right_top.second<<endl;
    clock_t start_time=clock();
    vector<POINT> Convex_Hull_point=func(point);
    clock_t end_time=clock();
    sort(point.begin(),point.end(),cmp_y_descend);
    //POINT right_bottom=point[n-1];
    //if (find(Convex_Hull_point.begin(),Convex_Hull_point.end(),right_bottom)==Convex_Hull_point.end()){
     //   Convex_Hull_point.push_back(right_bottom);
     //   Convex_Hull_point=GrahamScan(Convex_Hull_point);
    //}
    double duration=(double)(end_time-start_time)/(double)(CLOCKS_PER_SEC);
    string::size_type idx=file.find('.');
    cout << "The time of  "<<file.substr(0,idx)<<"  is:" << duration  <<"s"<<endl;
    CH(Convex_Hull_point,file.c_str());
}

int main()
{
    ofstream out1("point.txt");
    int n,range=101;
    cout<<"Please input the number of vertices:"<<endl;
    cin>>n;
    srand((int)time(NULL));

    data_generalization(point,n,range);
    sort(point.begin(),point.end());
    for (int i=0;i<n;i++){
        out1<<"("<<point[i].first<<","<<point[i].second<<")"<<endl;
    }
    Solve_solution_CH(BruteForce,point,"Brute_Force_Convex_Hull.txt");
    Solve_solution_CH(DivideConvex,point,"Divide_Conquer_Convex_Hull.txt");
    Solve_solution_CH(GrahamScan,point,"Graham_Scan_Convex_Hull.txt");

    return 0;
}

double multi(POINT &p1,POINT &p2,POINT &p0)
{   //p1代入p2p0直线
    return ((p1.first - p0.first)*(p2.second-p0.second) - (p1.second - p0.second)*(p2.first - p0.first));
}

double Distance(POINT p1,POINT p2)            //两个点的欧式距离
{
    return sqrt((p1.first-p2.first)*(p1.first-p2.first) + (p1.second-p2.second)*(p1.second-p2.second));
}

bool cmp_descend(POINT a,POINT b){
    if (a.first>b.first) return true;
    else return false;
}

bool cmp_y_descend(POINT &a, POINT &b)
{
    if (a.second<b.second) return true;
    else return false;

}
