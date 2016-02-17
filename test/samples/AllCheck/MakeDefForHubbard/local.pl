#!/usr/local/bin/perl
  #print "start !! \n";
  &input;
  #input!!
  #@intra_U=(1..5);
  #$inp=(1..10);
  $orb_num=$tmp_orb;
  $L_x=$tmp_Lx;
  $L_y=$tmp_Ly;
  $L=$L_x*$L_y;
  $Ns=$L;
  $All_N=$Ns*$orb_num;
  printf("CHECK local.prl L_x=$L_x L_y=$L_y All_N=$All_N orb=$orb_num \n");

  $fname="zlocspn.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "NlocalSpin 0  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========i_0LocSpn_1IteElc ====== \n";
  printf FILE "========misawa======== \n";

  for($i=0;$i<$All_N;$i+=1){
    printf FILE ("%4d  %4d\n",$i,1);
 }
 close(FILE);
 printf "lolal.prl finish \n"; 

#subroutine
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
  }
  if($Lx_cnt==0 || $Ly_cnt==0||$orb_cnt==0 ){
    printf "FAITAL ERROR IN input.txt !!!!!!!!!!!!! \n";
  }
  #input FINISH
 }

