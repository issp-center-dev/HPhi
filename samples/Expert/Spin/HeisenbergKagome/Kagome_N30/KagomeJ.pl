#!/usr/bin/perl -w 
  #input!!
  &input;
  #input!!
  $orb_num=$tmp_orb;
  $L_x=$tmp_Lx;
  $L_y=$tmp_Ly;
  $L=$L_x*$L_y;
  $Ns=$L;
  $All_N=$Ns*$orb_num;
  printf("CHECK J_inter $All_N  $J L_x=$L_x L_y=$L_y orb=$orb_num \n");

  $tmp=$All_N;
  $fname="zlocspn.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "loc $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========HundCoupling ====== \n";
  printf FILE "========misawa======== \n";

  for($i=0;$i<$All_N;$i+=1){
    printf FILE ("%4d  %4d\n",$i,0);
  }
  close(FILE);
 
# inter site J
  $tmp=2*$All_N;
  $fname="zcoulombinter.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "NHund $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========HundCoupling ====== \n";
  printf FILE "========misawa======== \n";

  $J_eff = -0.25*$J;
  &output;
  close(FILE);
 
# inter site J
  $tmp=2*$All_N;
  $fname="zhund.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "NHund $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========HundCoupling ====== \n";
  printf FILE "========misawa======== \n";

  $J_eff = -0.5*$J;
  &output;
  close(FILE);
 
# inter site J
  $tmp=2*$All_N;
  $fname="zexchange.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "NHund $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========HundCoupling ====== \n";
  printf FILE "========misawa======== \n";

  $J_eff = 0.5*$J;
  &output;
  close(FILE);


 
 sub input{
  #input START 
  $Lx_cnt=0;
  $Ly_cnt=0;
  $orb_cnt=0;
  $file=sprintf("input.txt");
  open(INPUTFILE,$file);
  while($name=<INPUTFILE>){
    chomp $name;
    #DELETE EMPTY
    $_=$name; 
    s/^\s+//;
    $name=$_; 
    @tmp = split /\s+/, $name;
    #printf "$tmp[0] $tmp[1] \n";
    if($tmp[0] eq 'Lx'){
      #printf "AA $tmp[0] $tmp[1] \n";
      $tmp_Lx = $tmp[1];
      $Lx_cnt=1;
    } 
    if($tmp[0] eq 'Ly'){
      #printf "AA $tmp[0] $tmp[1] \n";
      $tmp_Ly = $tmp[1];
      $Ly_cnt=1;
    } 
    if($tmp[0] eq 'orb_num'){
      #printf "AA $tmp[0] $tmp[1] \n";
      $tmp_orb = $tmp[1];
      $orb_cnt=1;
    } 
    if($tmp[0] eq 'J'){
      $J = $tmp[1];
    } 
  }
  if($Lx_cnt==0 || $Ly_cnt==0||$orb_cnt==0){
    printf "FAITAL ERROR IN input.txt !!!!!!!!!!!!! \n";
  }
  #input FINISH
 }

sub output{
   printf FILE ("%4d  %4d %lf\n",0,2,$J_eff);
   printf FILE ("%4d  %4d %lf\n",0,4,$J_eff);
   printf FILE ("%4d  %4d %lf\n",1,2,$J_eff);
   printf FILE ("%4d  %4d %lf\n",1,3,$J_eff);
   printf FILE ("%4d  %4d %lf\n",2,3,$J_eff);
   printf FILE ("%4d  %4d %lf\n",2,4,$J_eff);
   printf FILE ("%4d  %4d %lf\n",3,5,$J_eff);
   printf FILE ("%4d  %4d %lf\n",3,6,$J_eff);
   printf FILE ("%4d  %4d %lf\n",4,5,$J_eff);
   printf FILE ("%4d  %4d %lf\n",4,7,$J_eff);
   printf FILE ("%4d  %4d %lf\n",5,6,$J_eff);
   printf FILE ("%4d  %4d %lf\n",5,7,$J_eff);
   printf FILE ("%4d  %4d %lf\n",6,8,$J_eff);
   printf FILE ("%4d  %4d %lf\n",6,9,$J_eff);
   printf FILE ("%4d  %4d %lf\n",7,8,$J_eff);
   printf FILE ("%4d  %4d %lf\n",7,10,$J_eff);
   printf FILE ("%4d  %4d %lf\n",8,9,$J_eff);
   printf FILE ("%4d  %4d %lf\n",8,10,$J_eff);
   printf FILE ("%4d  %4d %lf\n",9,11,$J_eff);
   printf FILE ("%4d  %4d %lf\n",9,12,$J_eff);
   printf FILE ("%4d  %4d %lf\n",10,11,$J_eff);
   printf FILE ("%4d  %4d %lf\n",10,13,$J_eff);
   printf FILE ("%4d  %4d %lf\n",11,12,$J_eff);
   printf FILE ("%4d  %4d %lf\n",11,13,$J_eff);
   printf FILE ("%4d  %4d %lf\n",12,14,$J_eff);
   printf FILE ("%4d  %4d %lf\n",12,15,$J_eff);
   printf FILE ("%4d  %4d %lf\n",13,14,$J_eff);
   printf FILE ("%4d  %4d %lf\n",13,16,$J_eff);
   printf FILE ("%4d  %4d %lf\n",14,15,$J_eff);
   printf FILE ("%4d  %4d %lf\n",14,16,$J_eff);
   printf FILE ("%4d  %4d %lf\n",15,17,$J_eff);
   printf FILE ("%4d  %4d %lf\n",15,19,$J_eff);
   printf FILE ("%4d  %4d %lf\n",16,17,$J_eff);
   printf FILE ("%4d  %4d %lf\n",16,18,$J_eff);
   printf FILE ("%4d  %4d %lf\n",17,18,$J_eff);
   printf FILE ("%4d  %4d %lf\n",17,19,$J_eff);
   printf FILE ("%4d  %4d %lf\n",18,20,$J_eff);
   printf FILE ("%4d  %4d %lf\n",18,21,$J_eff);
   printf FILE ("%4d  %4d %lf\n",19,20,$J_eff);
   printf FILE ("%4d  %4d %lf\n",19,22,$J_eff);
   printf FILE ("%4d  %4d %lf\n",20,21,$J_eff);
   printf FILE ("%4d  %4d %lf\n",20,22,$J_eff);
   printf FILE ("%4d  %4d %lf\n",21,23,$J_eff);
   printf FILE ("%4d  %4d %lf\n",21,24,$J_eff);
   printf FILE ("%4d  %4d %lf\n",22,23,$J_eff);
   printf FILE ("%4d  %4d %lf\n",22,25,$J_eff);
   printf FILE ("%4d  %4d %lf\n",23,24,$J_eff);
   printf FILE ("%4d  %4d %lf\n",23,25,$J_eff);
   printf FILE ("%4d  %4d %lf\n",24,26,$J_eff);
   printf FILE ("%4d  %4d %lf\n",24,27,$J_eff);
   printf FILE ("%4d  %4d %lf\n",25,26,$J_eff);
   printf FILE ("%4d  %4d %lf\n",25,28,$J_eff);
   printf FILE ("%4d  %4d %lf\n",26,27,$J_eff);
   printf FILE ("%4d  %4d %lf\n",26,28,$J_eff);
   printf FILE ("%4d  %4d %lf\n",27,0,$J_eff);
   printf FILE ("%4d  %4d %lf\n",27,29,$J_eff);
   printf FILE ("%4d  %4d %lf\n",28,1,$J_eff);
   printf FILE ("%4d  %4d %lf\n",28,29,$J_eff);
   printf FILE ("%4d  %4d %lf\n",29,1,$J_eff);
   printf FILE ("%4d  %4d %lf\n",29,0,$J_eff);


}
