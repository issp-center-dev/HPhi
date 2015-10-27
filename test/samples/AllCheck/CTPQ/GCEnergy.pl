#!/usr/local/bin/perl
   $file=sprintf("GCCLanczos/output/zvo_CG_energy.dat");
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
     $Lanczos_energy = $foo[1];
     $cnt+=1;
   }
   close(INPUTFILE);

   $file=sprintf("tmp_TPQ.dat");
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
     $TPQ_energy = $foo[2];
     $cnt+=1;
   }
   close(INPUTFILE);

   $file=sprintf("tmp_FullDiag.dat");
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
     $Full_energy = $foo[1];
     $cnt+=1;
   }
   close(INPUTFILE);

   printf("\n");
   printf("Lanczos:   $Lanczos_energy\n");
   printf("TPQ:       $TPQ_energy\n");
   printf("FullDiag:  $Full_energy\n");
   printf("\n");

