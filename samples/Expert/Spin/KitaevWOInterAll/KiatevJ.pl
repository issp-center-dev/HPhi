#!/usr/bin/perl -w 
  #input!!
  &input;
  #input!!
  $orb_num = $tmp_orb;
  $L_x     = $tmp_Lx;
  $L_y     = $tmp_Ly;
  $L       = $L_x*$L_y;
  $Ns      = $L;
  $All_N   = $Ns*$orb_num;
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
    printf FILE ("%4d  %4d\n",$i,1);
  }
  close(FILE);
 
# inter site J
  $tmp=1*$All_N/2;
  $fname="zcoulombinter.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "NHund $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========HundCoupling ====== \n";
  printf FILE "========misawa======== \n";

  $J_eff   = -0.25*$J;
  &output_z;
  close(FILE);
 
# inter site J
  $tmp=1*$All_N/2;
  $fname="zhund.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "NHund $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========HundCoupling ====== \n";
  printf FILE "========misawa======== \n";

  $J_eff   = -0.5*$J;
  &output_z;
  close(FILE);
 
# inter site J
  $tmp=2*$All_N/2;
  $fname="zexchange.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "NHund $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========HundCoupling ====== \n";
  printf FILE "========misawa======== \n";

  $J_eff = 0.25*$J;
  &output_x;
  &output_y;

  close(FILE);
 
# inter site J
  $tmp=2*$All_N/2;
  $fname="zpairlift.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "NHund $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========HundCoupling ====== \n";
  printf FILE "========misawa======== \n";

  $J_eff = 0.25*$J;
  &output_x;
  $J_eff = -0.25*$J;
  &output_y;
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

sub output_x{
  for($all_i=0;$all_i<$All_N;$all_i+=1){
    $orb_i   = $all_i%$orb_num;
    $site_i  = ($all_i-$orb_i)/$orb_num;
    $x_i     = $site_i%$L_x;
    $y_i     = ($site_i-$x_i)/$L_x;
    if($orb_i==0){
      $orb_j = 1;
      $x_j   = $x_i;
      $y_j   = ($y_i+1)%$L_y;
      $all_j = ($x_j+$y_j*$L_x)*$orb_num+$orb_j;
      printf FILE ("%4d  %4d %lf\n",$all_i,$all_j,$J_eff);
    }
  }
}

sub output_y{
  for($all_i=0;$all_i<$All_N;$all_i+=1){
    $orb_i   = $all_i%$orb_num;
    $site_i  = ($all_i-$orb_i)/$orb_num;
    $x_i     = $site_i%$L_x;
    $y_i     = ($site_i-$x_i)/$L_x;
    if($orb_i==1){
      $orb_j = 0;
      $x_j   = ($x_i+1)%$L_x;
      $y_j   = $y_i;
      $all_j = ($x_j+$y_j*$L_x)*$orb_num+$orb_j;
      printf FILE ("%4d  %4d %lf\n",$all_i,$all_j,$J_eff);
    }
  }
}

sub output_z{
  for($all_i=0;$all_i<$All_N;$all_i+=1){
    $orb_i   = $all_i%$orb_num;
    $site_i  = ($all_i-$orb_i)/$orb_num;
    $x_i     = $site_i%$L_x;
    $y_i     = ($site_i-$x_i)/$L_x;
    if($orb_i==0){
      $orb_j = 1;
      $x_j   = $x_i;
      $y_j   = $y_i;
      $all_j = ($x_j+$y_j*$L_x)*$orb_num+$orb_j;
      printf FILE ("%4d  %4d %lf\n",$all_i,$all_j,$J_eff);
    }
  }
}
