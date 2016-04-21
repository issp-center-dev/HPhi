#!/usr/local/bin/perl
#===============def input====================
  $cnt =0;
  $file=sprintf("zInterAll-notdiag.def");
  open(INPUTFILE,$file);
  while($name=<INPUTFILE>){
    chomp $name;
    #DELETE EMPTY
    $_=$name; 
    s/^\s+//;
    $name=$_; 
    #DELETE EMPTY FINISH
    @foo = split /\s+/, $name;
    #printf "$cnt $foo[0] $foo[1] $foo[2] $foo[3] $foo[4] $foo[5] $foo[6] $foo[7] \n";
    if($cnt>4){
      $tmp_cnt=$cnt-5;
      for($i=0;$i<10;$i++){
        $phys[$tmp_cnt][$i]=$foo[$i];
      }
    }
    $cnt+=1;
  }
  close(INPUTFILE);
  $cnt_max=$cnt-5;
#=================================
 for($cnt=0;$cnt<$cnt_max;$cnt++){
   if($phys[$cnt][0]==$phys[$cnt][2] && $phys[$cnt][4]==$phys[$cnt][6]){
     printf("%d %d %d %d:cnt=%d: %f \n",$phys[$cnt][0],$phys[$cnt][2],$phys[$cnt][4],$phys[$cnt][6],$cnt,$phys[$cnt][8]);
   }
 } 
