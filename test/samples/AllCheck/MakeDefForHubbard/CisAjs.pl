#!/usr/local/bin/perl
  #input!!
  &input;
  #input!!
  $orb_num=$tmp_orb;
  $U=$tmp_U;
  $L_x=$tmp_Lx;
  $L_y=$tmp_Ly;
  $L=$L_x*$L_y;
  $Ns=$L;
  $All_N=$Ns*$orb_num;
  printf("CHECK CisAjs U=$U L_x=$L_x L_y=$L_y orb=$orb_num \n");
  $tmp=2*$All_N**2;
 
  $fname="zcisajs.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "NCisAjs $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========Green functions  ====== \n";
  printf FILE "========misawa======== \n";

  $all_cnt=0;
  for($i=0;$i<$All_N;$i+=1){
    for($j=0;$j<$All_N;$j+=1){
      printf FILE (" %4d %4d %4d %4d\n",$i, 0, $j,0);
      $all_cnt+=1;
    }
  }
  for($i=0;$i<$All_N;$i+=1){
    for($j=0;$j<$All_N;$j+=1){
      printf FILE (" %4d %4d %4d %4d\n",$i,1,$j,1);
      $all_cnt+=1;
    }
  }
  if($all_cnt!=2*$All_N**2){
    printf("ERROR IN CisAjs.prl \n");
  }
  printf("CisAjs.prl finish \n");

 sub input{
  #input START 
  $Lx_cnt=0;
  $Ly_cnt=0;
  $orb_cnt=0;
  $U_cnt=0;
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
    if($tmp[0] eq 'U'){
      #printf "AA $tmp[0] $tmp[1] \n";
      $tmp_U = $tmp[1];
      $U_cnt=1;
    } 
  }
  if($Lx_cnt==0 || $Ly_cnt==0||$orb_cnt==0){
    printf "FAITAL ERROR IN input.txt !!!!!!!!!!!!! \n";
  }
  #input FINISH
 }
 
