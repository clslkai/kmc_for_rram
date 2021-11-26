#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

    #define length 20
    #define width 75
    #define exp 2.71828
    #define freq 1e12
    const float ron=1;
    const float roff=1000;
    float voltage,temperature;
    float energy[12]={0.625,0.45,0.55,0.4,0.65,0.55,0.45,0.55,0.65,0.5,0.5,0.5};
    float rate[length][width][12];    //改成三维数组
    const float kb=8.617e-5;  //波尔兹曼常数，单位eV/K
    int ion_number,adatom_number,steps,h_filament,event=0;
    float model_time,time_step,total_rate;
    int array[length][width];
    float e_potential[length+1][width];
    float resistance[width],e_field[width];
    int i,j,k,x,y,flag,test;
    const float lattice_const=7.673e-10;
	int jo,growth;

int main()
{
   //下面是函数声明
    void ionstats();    //统计并随机选中一个离子，输出其坐标
    void printarray();  //输出数组
    void move();        //原子的移动
    void rate_calculate();
    int event_select();
    void oxide();
    void reduction();
    void desorption();
    void diffusion();
    void c_adsorption();
    void c_desorption();
    void c_reduction();
    void a_reduction();
    int get_end();
    void potential_calculate();
    
    //下面是变量的初始化
    temperature=300;
    model_time=0;
    voltage=3;
    srand((unsigned)time(0));
    flag=0;
	growth=0;
    //下面是器件的初始化
    for(j=0;j<width;j++)
      array[0][j]=1;
    for(i=1;i<length;i++)
      for(j=0;j<width;j++)
        array[i][j]=0;
    
    for(i=0;i<length+1;i++)
      for(j=0;j<width;j++)
        e_potential[i][j]=0;

    do
    //下面可以设置循环次数
   //  for(flag=0;flag<5e4;flag++)
     {
      flag=flag+1;      
//     for(j=0;j<width;j++)
//      {
//         i=length-1;
//         h_filament=0;
//         while(array[i][j]==1)
//          {
//            h_filament=h_filament+1;
//            i=i-1;
//          }
//            resistance[j]=ron*h_filament+roff*(length-h_filament);
//        }
        potential_calculate();
	    ionstats();
        //for(test=0;test<width;test++)
        //printf("e_field[%d]=%f ",test,e_field[test]);
        //printf("\n");
         event=event_select();
         printf("event=%d,flag=%d,x=%d,y=%d\n",event,flag,i,j);
         switch(event)
          {
           case 0:oxide(i,j);         break;
           case 1:desorption(i,j);    break;
           case 2:diffusion(i,j);     break;
           case 3:move(i,j);          break;
           case 4:c_reduction(i,j);   break;
           case 5:c_reduction(i,j);   break;
           case 6:c_reduction(i,j);   break;
           case 7:diffusion(i,j);     break;
           case 8:a_reduction(i,j);   break;
         }
      printarray();    //下面判断细丝百分比 

 // for(i=0;i<width;i++)
//	  {
 //      if(array[length-growth*3-2][i]==1)
//		{
 //       growth=growth+1;
 //   	 printf("growth=%d\n",growth);
//		 scanf("%d",&jo);
//		}
 //}

      }
     while(get_end());
 
 return(1);
}

void ionstats()   //统计离子情况，计算其活动概率,赋值给rate矩阵
{
  //下面是概率矩阵的初始化
 for(i=0;i<length;i++)
   for(j=0;j<width;j++)
     for(k=0;k<12;k++)
      rate[i][j][k]=0;
    
  for(j=0;j<width;j++)
    {
      //氧化
      //下面讨论器件顶层离子
      if(array[1][j]==0)
       rate[1][j][0]=freq*pow(exp,(energy[0]-0.5*(e_potential[0][j]-e_potential[1][j]))/(kb*temperature*(-1)));
        
      if(array[1][j]==2)
        {
          switch(array[2][j])
          {
          case 2:
               {
                  rate[1][j][2]=freq*pow(exp,energy[2]/(kb*temperature*(-1)));
                  rate[1][j][8]=freq*pow(exp,(energy[8]+0.5*(e_potential[0][j]-e_potential[1][j]))/(kb*temperature*(-1)));
                  break;
               } 
          case 1:
               {
                  rate[1][j][2]=freq*pow(exp,energy[2]/(kb*temperature*(-1)));
                  rate[1][j][4]=freq*pow(exp,(energy[4]-0.5*(e_potential[0][j]-e_potential[1][j]))/(kb*temperature*(-1)));
                  rate[1][j][8]=freq*pow(exp,(energy[8]+0.5*(e_potential[0][j]-e_potential[1][j]))/(kb*temperature*(-1)));
                  break;
               }
          case 0:
               {
                  rate[1][j][1]=freq*pow(exp,(energy[1]-0.5*(e_potential[1][j]-e_potential[2][j]))/(kb*temperature*(-1)));
                  rate[1][j][2]=freq*pow(exp,energy[2]/(kb*temperature*(-1)));
                  rate[1][j][8]=freq*pow(exp,(energy[8]+0.5*(e_potential[0][j]-e_potential[1][j]))/(kb*temperature*(-1)));
                  break;
               }
          }       
       }
    }
  for(i=2;i<length;i++)
   for(j=0;j<width;j++)
   {
     if(array[i][j]==2)  //器件内离子
      {
        if(i==length-1||array[i+1][j]==1||array[i][j-1]==1||array[i][j+1]==1)   //可能被还原
           {
           if(array[i][j-1]==1&&array[i][j+1]==1&&(i==length-1||array[i+1][j]==1))   //两侧都是铜原子
               rate[i][j][6]=freq*pow(exp,(energy[6]-0.5*(e_potential[i-1][j]-e_potential[i][j]))/(kb*temperature*(-1)));
           else if((array[i][j-1]==1&&(i==length-1||array[i+1][j]==1))||(array[i][j+1]==1&&(i==length-1||array[i+1][j]==1)))
                 rate[i][j][5]=freq*pow(exp,(energy[5]-0.5*(e_potential[i-1][j]-e_potential[i][j]))/(kb*temperature*(-1)));
                else
                 {
                  rate[i][j][4]=freq*pow(exp,(energy[4]-0.5*(e_potential[i-1][j]-e_potential[i][j]))/(kb*temperature*(-1)));
                  rate[i][j][7]=freq*pow(exp,energy[7]/(kb*temperature*(-1)));
                 }
           }
          else
           rate[i][j][3]=freq*pow(exp,energy[3]/(kb*temperature*(-1)));
      }
   }
}

int event_select()
{
  float event_pick;
  int h;
  total_rate=0;
 for(i=0;i<length;i++)
    for(j=0;j<width;j++)
       for(h=0;h<12;h++)
          {
           if(rate[i][j][h]!=0)
             {
           //  printf("rate[%d][%d][%d]=%f ",i,j,h,rate[i][j][h]); 
               total_rate=total_rate+rate[i][j][h];
             }
          }
  event_pick=((float)rand()/RAND_MAX)*total_rate;
  time_step=1/total_rate*1e9;
  model_time=model_time+time_step/1e3;
  printf("Total Rate=%f   ",total_rate);
  printf("Event Pick=%f  ",event_pick);
  printf("Time Step=%fns  ",time_step);
  printf("Switching Time=%fus\n",model_time);
  for(i=0;i<length;i++)
     for(j=0;j<width;j++)
          for(h=0;h<12;h++)
          {
      // printf("rate[%d][%d][%d]=%f ",i,j,h,rate[i][j][h]); 
         if(event_pick<rate[i][j][h])
           return h;
          else 
           event_pick=event_pick-rate[i][j][h];
          }
}


void printarray()/*输出数组*/
{
int i,j;
for(i=0;i<length;i++)/*输出数组的所有的值*/
{
for(j=0;j<width;j++)
printf("%d",array[i][j]);
printf("\n");
}
}


void move(int a,int b) //离子电介质内移动
{
int h,path;
float dir2[4];        //每个区间对应的一个数值
float dir;            //概率的发生在哪个区间
// printf("a=%d,b=%d\n",a,b);
for(h=0;h<4;h++)
 dir2[h]=0;
 if(array[a-1][b]==0)
 dir2[0]=freq*pow(exp,(energy[5]+0.5*(e_potential[a-1][b]-e_potential[a][b]))/(kb*temperature*(-1)));            //上
 if(a!=(length-1)&&array[a+1][b]==0)
 dir2[1]=dir2[0]+freq*pow(exp,(energy[5]+0.5*(e_potential[a+1][b]-e_potential[a][b]))/(kb*temperature*(-1)));   //下
 else dir2[1]=dir2[0];
 if(array[a][b-1]==0)
 dir2[2]=dir2[1]+freq*pow(exp,(energy[5]+0.5*(e_potential[a][b-1]-e_potential[a][b]))/(kb*temperature*(-1)));       //左
 else dir2[2]=dir2[1];
 if(array[a][b+1]==0)
 dir2[3]=dir2[2]+freq*pow(exp,(energy[5]+0.5*(e_potential[a][b+1]-e_potential[a][b]))/(kb*temperature*(-1)));     //右
 else  dir2[3]=dir2[2];
array[a][b]=0;
dir=(float)rand()/RAND_MAX*dir2[3];
if (dir<=dir2[0])
path=1;
else if (dir>dir2[0]&&dir<=dir2[1])
path=2;
else if (dir>dir2[1]&&dir<=dir2[2])
path=3;
else if (dir>dir2[2]&&dir<=dir2[3])
path=4;
else
printf("error\n");
printf("path=%d\n",path);

switch (path)
{
case 1:a=a-1;          break;
case 2:a=a+1;          break;
case 3:b=(b+1)%width;  break;
case 4:b=(b-1+width)%width;  break;
default:printf("error\n");

}
array[a][b]=2;
//printf("a=%d,b=%d",a,b);

}


void oxide(int a, int b)
{
array[a][b]=2;
}

void reduction(int a,int b)
{
array[a][b]=0;
}

void desorption(int a, int b)
{
array[a][b]=0;
array[a+1][b]=2;
}

void diffusion(int a,int b)
{
int h,path;
float dir2[2];        //每个区间对应的一个数值
float dir;            //概率的发生在哪个区间
// printf("a=%d,b=%d\n",a,b);
for(h=0;h<2;h++)
 dir2[h]=0;
 dir2[0]=dir2[0]+freq*pow(exp,(energy[2]+0.5*(e_potential[a][b-1]-e_potential[a][b]))/(kb*temperature*(-1)));       //左
 dir2[1]=dir2[1]+freq*pow(exp,(energy[2]+0.5*(e_potential[a][b+1]-e_potential[a][b]))/(kb*temperature*(-1)));     //右
array[a][b]=0;
dir=(float)rand()/RAND_MAX*dir2[1];
if (dir<=dir2[0])
path=1;
else if (dir>dir2[0]&&dir<=dir2[1])
path=2;
else
printf("error\n");
//printf("path=%d\n",path);

switch (path)
{
case 1:b=(b+1)%width;  break;
case 2:b=(b-1+width)%width;  break;
default:printf("error\n");

}
array[a][b]=2;
//printf("a=%d,b=%d",a,b);

}

void c_adsorption(int a,int b)
{
array[a][b]=0;
array[a+1][b]=2;
}
void c_desorption(int a,int b)
{
array[a][b]=0;
array[a-1][b]=2;
}
void c_reduction(int a,int b)
{
array[a][b]=1;
}

void a_reduction(int a,int b)
{
array[a][b]=0;
}

int get_end()
{

  int a;
  int k;
  k=1;
  for(a=0;a<width;a++)
  {
   if(array[1][a]==1)
   k=0; 
  }
return(k);
} 

void potential_calculate()
{
  //用有限差分法求解矩形域上的Poisson方程
    float z,tol,norm;
    tol=1e-4;
    norm=0;
    z=0;
    for(j=0;j<width;j++)
      {
      e_potential[0][j]=voltage;
      e_potential[length][j]=0;
      }
    for(i=1;i<length;i++)
      for(j=0;j<width;j++)
       {
         if(array[i][j]==1)     
           e_potential[i][j]=0;
       }
    do    
    {
       norm=0;
       for(i=1;i<length;i++)
         for(j=0;j<width;j++)
           {
             if(array[i][j]!=1)
               {
                 if(j==0)
                 z=(e_potential[i-1][j]+e_potential[i+1][j]+e_potential[i][j+1])/3;
                 else if(j==(width-1))
                 z=(e_potential[i-1][j]+e_potential[i+1][j]+e_potential[i][j-1])/3;
                 else  
                 z=(e_potential[i-1][j]+e_potential[i+1][j]+e_potential[i][j-1]+e_potential[i][j+1])/4;
                 if(fabs(e_potential[i][j]-z)>norm)
                    {
                     norm=fabs(e_potential[i][j]-z);
                     e_potential[i][j]=z;
                    }
               }  
           }
    } 
   while(norm>tol);
 // for(i=1;i<length;i++)/*输出数组的所有的值*/
 //  {
 //    for(j=1;j<width+1;j++)
 //       printf("%3.2f ",e_potential[i][j]);
 //       printf("\n");
 //  }
}
