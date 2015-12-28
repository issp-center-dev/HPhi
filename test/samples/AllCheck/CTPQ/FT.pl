#!/usr/local/bin/perl

#===============def input====================
    $min_ene=0;
    $cnt=0;
    $file=sprintf("Eigenvalue.dat ");
    open(INPUTFILE,$file);
    while($name=<INPUTFILE>){
      chomp $name;
      #DELETE EMPTY
      $_=$name; 
      s/^\s+//;
      $name=$_; 
      #DELETE EMPTY FINISH
      @foo = split /\s+/, $name;
      #printf "$cnt $foo[0] $foo[1] $foo[2] $foo[3] $foo[4] \n";
      $Ene[$cnt] = $foo[1];  
      $cnt+=1;
    }
    close(INPUTFILE);
    $cnt_max=$cnt;
#================Sq=================
    $fname=sprintf("FT_phys.dat");
    open(FILE,">$fname");
    printf FILE "# beta E E^2\n";
    $delta_beta =0.01;  
    for($beta=0.01;$beta<100;$beta+=$delta_beta){
      if($beta>=0.1){
        $delta_beta =0.1;  
      }elsif($beta>=1){
        $delta_beta =1;  
      }
#      printf("$beta \n");
      $Z  = 0.0;
      $E  = 0.0;
      $E2 = 0.0;
      for($cnt=0;$cnt<$cnt_max;$cnt++){
        $tmp = $Ene[$cnt]-$min_ene;
        $Z  += exp(-$beta*$tmp);
        $E  += $Ene[$cnt]*exp(-$beta*$tmp);
        $E2 += $Ene[$cnt]*$Ene[$cnt]*exp(-$beta*$tmp);
      }
      $E  = $E/$Z;
      $E2 = $E2/$Z;
      $temperature = 1.0/$beta;
      printf FILE "$temperature $E $E2 \n";
    } 
    close(FILE);
