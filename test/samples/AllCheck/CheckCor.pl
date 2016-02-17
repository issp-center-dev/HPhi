#!/usr/local/bin/perl
   $L_x = 6;
   $file=sprintf("%s",$ARGV[0]);
   #print "$file \n";
   $cnt = 0;
   open(INPUTFILE,$file);
   while($name=<INPUTFILE>){
     chomp $name;
     #DELETE EMPTY
     $_=$name; 
     s/^\s+//;
     $name=$_; 
     #DELETE EMPTY FINISH
     @foo = split /\s+/, $name;
     $site_1 = $foo[0];
     $site_2 = $foo[1];
     $spn_1  = $foo[2];
     $site_3 = $foo[3];
     $site_4 = $foo[4];
     $spn_2  = $foo[5];
     $Cor[$site_1][$site_2][$spn_1][$site_3][$site_4][$spn_2]   = $foo[6];
     $cnt+=1;
   }
   close(INPUTFILE);


#==== nSS ====
   $nn_S  = 0.0;
   $nnn_S = 0.0;
   $site  = 0;
#========================
   $x     = $site%$L_x;
# nn
   $nx_1   =($x+1+$L_x)%$L_x;
   $nsite  = $nx_1;
   &add_nn_S;
   printf("nn spin cor. = $nn_S \n");
# nnn
   $nx_1   =($x+2+$L_x)%$L_x;
   $nsite  = $nx_1;
   &add_nnn_S;
   printf("nnn spin cor. = $nnn_S \n");
 sub add_nn_S{
    $nn_S  += $Cor[$site][$site][0][$nsite][$nsite][0]/4.0;
    $nn_S  += -$Cor[$site][$site][0][$nsite][$nsite][1]/4.0;
    $nn_S  += -$Cor[$site][$site][1][$nsite][$nsite][0]/4.0;
    $nn_S  += $Cor[$site][$site][1][$nsite][$nsite][1]/4.0;
    $nn_S  += -$Cor[$site][$nsite][0][$nsite][$site][1]/2.0;
    $nn_S  += -$Cor[$site][$nsite][1][$nsite][$site][0]/2.0;
 }

 sub add_nnn_S{
    $nnn_S  += $Cor[$site][$site][0][$nsite][$nsite][0]/4.0;
    $nnn_S  += -$Cor[$site][$site][0][$nsite][$nsite][1]/4.0;
    $nnn_S  += -$Cor[$site][$site][1][$nsite][$nsite][0]/4.0;
    $nnn_S  += $Cor[$site][$site][1][$nsite][$nsite][1]/4.0;
    $nnn_S  += -$Cor[$site][$nsite][0][$nsite][$site][1]/2.0;
    $nnn_S  += -$Cor[$site][$nsite][1][$nsite][$site][0]/2.0;
 }
