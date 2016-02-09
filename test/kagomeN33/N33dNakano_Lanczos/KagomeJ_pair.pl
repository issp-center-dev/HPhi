#!/usr/bin/perl -w 
  #input!!
  &input;
  &pairinput;
  #input!!
  $orb_num=$tmp_orb;
  $L_x=$tmp_Lx;
  $L_y=$tmp_Ly;
  $L=$L_x*$L_y;
  $Ns=$L;
  $All_N=$Ns*$orb_num;
  printf("CHECK J_inter $All_N  $J L_x=$L_x L_y=$L_y orb=$orb_num \n");

  $tmp=$All_N*2;
  $fname="greenone.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "CisAjs $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========mag ====== \n";
  printf FILE "========misawa======== \n";
  for($i=0;$i<$All_N;$i+=1){
    printf FILE ("%4d  %4d %4d  %4d\n",$i,0,$i,0);
    printf FILE ("%4d  %4d %4d  %4d\n",$i,1,$i,1);
  }
  close(FILE);


  $tmp=$All_N;
  $fname="zlocspn.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "loc $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========HundCoupling ====== \n";
  printf FILE "========misawa======== \n";

  for($i=0;$i<$All_N;$i+=1){
    printf FILE ("%4d  %4d\n",$i,1);
  }
  close(FILE);

#magnetic field
  $tmp=$All_N*2;
  $fname="zTrans.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "loc $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========HundCoupling ====== \n";
  printf FILE "========misawa======== \n";
  for($all_i=0;$all_i<$All_N;$all_i+=1){
    # 0.5 -> H*Sz
    printf FILE ("%4d %4d %4d %4d %lf 0.0\n",$all_i,0,$all_i,0,0.5*$H);
    printf FILE ("%4d %4d %4d %4d %lf 0.0\n",$all_i,1,$all_i,1,-0.5*$H);
  }
  close(FILE);
 
# inter site J
  $tmp=(4+2+2)*$cnt_max;
  $fname="zInterAll.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "NHund $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========HundCoupling ====== \n";
  printf FILE "========misawa======== \n";

  $J_eff = -0.25*$J; # coulomb inter
  for($cnt=0;$cnt<$cnt_max;$cnt++){
    for($spn1=0;$spn1<2;$spn1++){
      for($spn2=0;$spn2<2;$spn2++){
        printf FILE ("%4d %4d %4d %4d %4d %4d %4d %4d %lf 0.0\n",$siteI[$cnt],$spn1,$siteI[$cnt],$spn1,$siteJ[$cnt],$spn2,$siteJ[$cnt],$spn2,$J_eff);
      }
    } 
  }
  #&output;
  $J_eff = 0.5*$J;  # Hund
  for($cnt=0;$cnt<$cnt_max;$cnt++){
    for($spn1=0;$spn1<2;$spn1++){
      printf FILE ("%4d %4d %4d %4d %4d %4d %4d %4d %lf 0.0\n",$siteI[$cnt],$spn1,$siteI[$cnt],$spn1,$siteJ[$cnt],$spn1,$siteJ[$cnt],$spn1,$J_eff);
    } 
  }
  #&output;
  $J_eff =  0.5*$J;   # exchange
  for($cnt=0;$cnt<$cnt_max;$cnt++){
    for($spn1=0;$spn1<2;$spn1++){
      $spn2=1-$spn1; 
      printf FILE ("%4d %4d %4d %4d %4d %4d %4d %4d %lf 0.0\n",$siteI[$cnt],$spn1,$siteI[$cnt],$spn2,$siteJ[$cnt],$spn2,$siteJ[$cnt],$spn1,$J_eff);
    } 
  }
 
  #&output;
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
    if($tmp[0] eq 'H'){
      $H = $tmp[1];
    } 
  }
  if($Lx_cnt==0 || $Ly_cnt==0||$orb_cnt==0){
    printf "FAITAL ERROR IN input.txt !!!!!!!!!!!!! \n";
  }
  #input FINISH
 }
 
 
 sub pairinput{
  #input START 
  $cnt=0;
  $file=sprintf("pairlist.txt");
  open(INPUTFILE,$file);
  while($name=<INPUTFILE>){
    chomp $name;
    #DELETE EMPTY
    $_=$name; 
    s/^\s+//;
    $name=$_; 
    @tmp = split /\s+/, $name;
    $siteI[$cnt]=$tmp[0];
    $siteJ[$cnt]=$tmp[1];
    printf "$cnt $tmp[0] $tmp[1] \n";
    $cnt+=1;
  }
  $cnt_max=$cnt;
  #input FINISH
 }
