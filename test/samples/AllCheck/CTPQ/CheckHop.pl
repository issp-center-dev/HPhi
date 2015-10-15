#!/usr/local/bin/perl
   $cnt=0;
   $file=sprintf("%s",$ARGV[0]);
   open(INPUTFILE,$file);
   while($name=<INPUTFILE>){
     chomp $name;
     #DELETE EMPTY
     $_=$name; 
     s/^\s+//;
     $name=$_; 
     #DELETE EMPTY FINISH
     @foo = split /\s+/, $name;
     #printf "$cnt $foo[0] $foo[1] $foo[2] $foo[3] $foo[4] $foo[5] $foo[6]  \n";
     #for($i=0;$i<7;$i++){
     #  $phys[$cnt][$i]=$foo[$i];
     #}
     if($foo[0]==0 && $foo[2]==1){
       print("spin $foo[1]: hopping 0 1 = $foo[4] \n"); 
     }
     if($foo[0]==1 && $foo[2]==0){
       print("spin $foo[1]: hopping 1 0 = $foo[4] \n"); 
     }
     if($foo[0]==0 && $foo[2]==3){
       print("spin $foo[1]: hopping 0 3 = $foo[4] \n"); 
     }
     if($foo[0]==3 && $foo[2]==0){
       print("spin $foo[1]: hopping 3 0 = $foo[4] \n"); 
     }
     $cnt+=1;
   }
   close(INPUTFILE);

