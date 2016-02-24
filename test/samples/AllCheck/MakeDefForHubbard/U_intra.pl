#!/usr/local/bin/perl
  #print "start !! \n";
  #input!!
  &input;
  #input!!
  #@intra_U=(1..5);
  #$inp=(1..10);
  $U=$tmp_U;
  $orb_num=$tmp_orb;
  $L_x=$tmp_Lx;
  $L_y=$tmp_Ly;
  $L=$L_x*$L_y;
  $Ns=$L;
  $All_N=$Ns*$orb_num;
  printf("CHECK U_intra L_x=$L_x L_y=$L_y  All_N=$All_N orb=$orb_num \n");
 
  $tmp=$Ns*$orb_num;
  $fname="zcoulombintra.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "NCoulombIntra $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========i_0LocSpn_1IteElc ====== \n";
  printf FILE "========misawa======== \n";

  for($site_i=0;$site_i<$Ns;$site_i+=1){
    for($orb_i=0;$orb_i<$orb_num;$orb_i++){
      $all_i=$orb_num*$site_i+$orb_i;
      printf FILE ("%4d  %lf\n",$all_i,$U);
    }
  }

  close(FILE);
  printf "U_intra.prl finish \n";
#subroutine
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
  if($Lx_cnt==0 || $Ly_cnt==0||$orb_cnt==0 || $U_cnt==0){
    printf "FAITAL ERROR IN input.txt !!!!!!!!!!!!! \n";
  }
  #input FINISH
 }


