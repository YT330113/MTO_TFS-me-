Info<<"cost function"<<endl;
gradD=fvc::grad(D);  
sigmaD=mu*twoSymm(gradD) + (lambda*I)*tr(gradD);  
volScalarField Energy(0.25*(gradD+gradD.T())&&sigmaD); 

C=fvc::domainIntegrate(Energy).value()/C0;
MeanT=fvc::domainIntegrate(T).value()/area;
V=fvc::domainIntegrate(xh).value()/area-voluse;  

PowerDiss=0;
for (i=0; i<nObjPatch; i++)
{
     PowerDiss =PowerDiss - sum(
     phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
     + 0.5*magSqr(U.boundaryField()[conPatchList[i]]))); 
}
PowerDiss=PowerDiss/PowerDiss0;
if(GeoDim==2)
{
 PowerDiss=PowerDiss/std::pow(mesh.V()[0],1.0/3);
}
reduce(PowerDiss, sumOp<scalar>()); 

//输出目标及约束值
if (Pstream::master())
{ 
   ofstream outfile1("VoluseChange.txt",std::ios::app);           // 创建一个ofstream对象outfile1来管理输出流，app表示文件模式为追加到文件尾
   outfile1<<runTime.timeName()<<"    "<<V<<"\n";
   outfile1.close();

   ofstream outfile2("PowerDissRatio.txt",std::ios::app);
   outfile2<<runTime.timeName()<<"    "<<PowerDiss<<"\n";
   outfile2.close(); 

   ofstream outfile3("MeanT.txt",std::ios::app);
   outfile3<<runTime.timeName()<<"    "<<MeanT<<"\n";
   outfile3.close(); 

   ofstream outfile4("Time.txt",std::ios::app);
   outfile4<<runTime.elapsedClockTime()<<"\n";
   outfile4.close();   

   ofstream outfile5("ComplianceRatio.txt",std::ios::app);
   outfile5<<runTime.timeName()<<"    "<<C<<"\n";
   outfile5.close();  
}
