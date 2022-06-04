//Author: Yu Minghao    Updated: 2021.8.28 

static char help[] = "topology optimization of thermal-fluid-structural problem\n";
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"//用于粘度模型的单向输运模型，用于不可压单向流
#include "turbulentTransportModel.H"
#include "simpleControl.H"//simple算法的相关控制，用于计算过程中的收敛性判断
#include "fvOptions.H"//
#include <MMA.h>
#include <diff.c>

int main(int argc, char *argv[])
{
    PetscInitialize(&argc,&argv,PETSC_NULL,help);
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFvOptions.H"
    #include "createFields.H"//一些物理场的声明
    #include "readTransportProperties.H" //一些物理场的声明
    #include "readThermalProperties.H" //一些物理场的声明
    #include "readMechanicalProperties.H" //一些物理场的声明
     #include "initContinuityErrs.H"// 用于声明与初始化连续性误差，初始误差为0
    #include "opt_initialization.H"//一些参数的声明和初始化
    while (simple.loop(runTime))
    {
        #include "update.H"//根据设计变量更新材料属性
        #include "NS.H"//NS方程求解
        #include "HeatTransfer.H"//传热方程求解
        #include "LinearElasticity.H"//平衡方程求解
        #include "AdjHeatTransfer.H"//伴随能量方程(of 平均温度)
        #include "AdjNS_HT.H"//伴随NS方程(of 平均温度)
        #include "AdjNS_FF.H"//伴随能量方程(of 能量耗散)
        #include "costfunction.H"//计算目标和约束函数值             
        #include "sensitivity.H"//灵敏度分析及MMA
    }
    #include "finalize.H"//程序终止，delete变量
    return 0;
}
