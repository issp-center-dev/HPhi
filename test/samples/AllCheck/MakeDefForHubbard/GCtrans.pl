#!/usr/local/bin/perl
  #print "start !! \n";
  #input!!
  &input;
  #input!!
  $orb_num=$tmp_orb;
  $L_x=$tmp_Lx;
  $L_y=$tmp_Ly;
  $L=$L_x*$L_y;
  $Ns=$L;
  $All_N=$Ns*$orb_num;
  printf("CHECK trans $U L_x=$L_x L_y=$L_y  orb=$orb_num \n");
  $trans_1 = -1.0;
  $trans_2 = $t2;

  $tmp=3*$orb_num*$All_N;
  $tmp=2*$tmp;
  $fname="GCztransfer.def";
  open(FILE,">$fname");
  printf FILE "========misawa======== \n";
  printf FILE "NTransfer $tmp  \n";
  printf FILE "========misawa======== \n";
  printf FILE "========i_j_s_tijs  ====== \n";
  printf FILE "========misawa======== \n";

  # up   spin_i = 0
  # down spin_i = 1
for($spin_i=0;$spin_i<=1;$spin_i++){
  for($all_i=0;$all_i<$All_N;$all_i+=1){
    $orb_i  = $all_i%$orb_num;
    $site_i = ($all_i-$orb_i)/$orb_num;
    $i_x    = $site_i%$L_x;
    $i_y    = ($site_i-$i_x)/$L_x;

    #printf FILE ("%4d  %4d  %4d  %f \n",$all_i,$all_i,$spin_i,$U/2);
    printf FILE ("%4d %4d  %4d  %4d  %f %f\n",$all_i,$spin_i,$all_i,$spin_i,$U/2,0.0);

    #nearest
    #(1,0)
    $tmp_x    = ($i_x+1)%$L_x;
    $tmp_y    = $i_y;
    $tmp_site = $tmp_y*$L_x+$tmp_x;
    $orb_j    = 0;
    $all_j    = $tmp_site*$orb_num+$orb_j;
    $tmp_trans= $trans_1;
    #if($i_x==$L_x-1){
    #  printf FILE ("%4d  %4d  %4d  %f \n",$all_i,$all_j,$spin_i,$tmp_trans);
    #}else{
      printf FILE ("%4d  %4d  %4d  %4d %f %f\n",$all_i,$spin_i,$all_j,$spin_i,-$tmp_trans,0);
    #}
    #(-1,0)
    $tmp_x    = ($i_x-1+$L_x)%$L_x;
    $tmp_y    = $i_y;
    $tmp_site = $tmp_y*$L_x+$tmp_x;
    $orb_j    = 0;
    $all_j    = $tmp_site*$orb_num+$orb_j;
    $tmp_trans= $trans_1;
    #if($i_x==0){
    #  printf FILE ("%4d  %4d  %4d  %f \n",$all_i,$all_j,$spin_i,$tmp_trans);
    #}else{
      printf FILE ("%4d  %4d  %4d  %4d %f %f\n",$all_i,$spin_i,$all_j,$spin_i,-$tmp_trans,0);
    #}

    #(0,1)
    $tmp_x    = ($i_x+0)%$L_x;
    $tmp_y    = ($i_y+1)%$L_y;
    $tmp_site = $tmp_y*$L_x+$tmp_x;
    $orb_j    = 0;
    $all_j    = $tmp_site*$orb_num+$orb_j;
    $tmp_trans= $trans_1;
    #printf FILE ("%4d  %4d  %4d  %f \n",$all_i,$all_j,$spin_i,-$tmp_trans);

    #(0,-1)
    $tmp_x    = ($i_x+0)%$L_x;
    $tmp_y    = ($i_y-1+$L_y)%$L_y;
    $tmp_site = $tmp_y*$L_x+$tmp_x;
    $orb_j    = 0;
    $all_j    = $tmp_site*$orb_num+$orb_j;
    $tmp_trans= $trans_1;
    #printf FILE ("%4d  %4d  %4d  %f \n",$all_i,$all_j,$spin_i,-$tmp_trans);

    #next nearest
    #(1,1)
    $tmp_x    = ($i_x+1)%$L_x;
    $tmp_y    = ($i_y+1)%$L_y;
    $tmp_site = $tmp_y*$L_x+$tmp_x;
    $orb_j    = 0;
    $all_j    = $tmp_site*$orb_num+$orb_j;
    $tmp_trans=$trans_2;

    #printf FILE ("%4d  %4d  %4d  %f\n",$all_i,$all_j,$spin_i,-$tmp_trans);

    #(-1,-1)
    $tmp_x    = ($i_x-1+$L_x)%$L_x;
    $tmp_y    = ($i_y-1+$L_y)%$L_y;
    $tmp_site = $tmp_y*$L_x+$tmp_x;
    $orb_j    = 0;
    $all_j    = $tmp_site*$orb_num+$orb_j;
    $tmp_trans= $trans_2;

    #printf FILE ("%4d  %4d  %4d  %f\n",$all_i,$all_j,$spin_i,-$tmp_trans);

    #(1,-1)
    $tmp_x    = ($i_x+1)%$L_x;
    $tmp_y    = ($i_y-1+$L_y)%$L_y;
    $tmp_site = $tmp_y*$L_x+$tmp_x;
    $orb_j    = 0;
    $all_j    = $tmp_site*$orb_num+$orb_j;
    $tmp_trans=$trans_2;

    #printf FILE ("%4d  %4d  %4d  %f\n",$all_i,$all_j,$spin_i,-$tmp_trans);

    #(-1,1)
    $tmp_x    = ($i_x-1+$L_x)%$L_x;
    $tmp_y    = ($i_y+1+$L_y)%$L_y;
    $tmp_site = $tmp_y*$L_x+$tmp_x;
    $orb_j    = 0;
    $all_j    = $tmp_site*$orb_num+$orb_j;
    $tmp_trans= $trans_2;

    #printf FILE ("%4d  %4d  %4d  %f\n",$all_i,$all_j,$spin_i,-$tmp_trans);
 
  }
}
 close(FILE);
 printf "trans.prl finish \n";

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
    if($tmp[0] eq 'U'){
      #printf "AA $tmp[0] $tmp[1] \n";
      $U  = $tmp[1];
    } 
    if($tmp[0] eq 't2'){
      #printf "AA $tmp[0] $tmp[1] \n";
      $t2  = $tmp[1];
    } 
  }
  if($Lx_cnt==0 || $Ly_cnt==0||$orb_cnt==0){
    printf "FAITAL ERROR IN input.txt !!!!!!!!!!!!! \n";
  }
  #input FINISH
 }


