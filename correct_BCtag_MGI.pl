#!/usr/bin/env perl
use warnings;

if (@ARGV != 5) {  
  print
"\n#####################################################################################
Usage: perl $0 <inbam> <outbam> <BCbinmap1> <BCbinmap2> <samtools-executable>  \n
Please drop your suggestions and clarifications to <christoph.ziegenhain\@ki.se>\n
######################################################################################\n\n";
  exit;
}
BEGIN{
$inbam = $ARGV[0];  
$outbam = $ARGV[1];  
$binmap1 = $ARGV[2];  
$binmap2 = $ARGV[3];  
$samtoolsexc = $ARGV[4]; 
}
$argLine = join(" ", @ARGV);

open(INBAM, "$samtoolsexc view -h $inbam | sed 's/BC:Z://' | " ) || die "Couldn't open file $inbam. Check permissions!\n Check if it is a bam file and it exists\n\n";

open(BCBAM,"| $samtoolsexc view -b -o $outbam -");
$bamhead = 0;

my %bcmap1; 
open(DATA1, "cat $binmap1 | sed 's/,/\t/g' | cut -f1,3 | grep -v 'falseBC' | ") || die "Can't open $binmap1 ! \n";
while (<DATA1>) {
  my ($raw, @fixedBC) = split(/\t/);
  $bcmap1{$raw} = \@fixedBC;
}
close DATA1;

my %bcmap2;  
open(DATA2, "cat $binmap2 | awk -F',' '\$2==0 {print \$1,\$3}' | ") || die "Can't open $binmap2!\n";  
while (<DATA2>) {  
    my ($falseBC, $trueBC) = split(/\s+/);  
    $bcmap2{$falseBC} = $trueBC;  
}  
close DATA2; 


while (<INBAM>) {
  $read=$_;

 if($read =~ /^\@/){
   print BCBAM $read;
   next;
 }

 chomp($read);
  @read = split(/\t/,$read);
  $thisBC = $read[11];

  if (defined($bcmap1{$thisBC})) {
    #print "BC is in hash\n";
    $correctBC = $bcmap1{$thisBC}[0];
    chomp($correctBC);
  }
  else {
    #print "BC is not in hash\n";
    $correctBC = $thisBC;
  }

  if (defined($bcmap2{$correctBC})) {  
        $correctBCtrans = $bcmap2{$correctBC};  
  } else {  
        $correctBCtrans = $correctBC;  
  } 

        if(!$bamhead){
          print(BCBAM join("\t", ("@"."PG","ID:zUMIs-fqfilter","PN:zUMIs-correct_BCtag", "VN:2","CL:correct_BCtag.pl ${argLine}")) . "\n");
          $bamhead = 1;
        }
  #print STDERR "Debug: correctBCtrans = $correctBCtrans\n";

  if ($read[1] == 77 && $correctBCtrans eq $correctBC) {  
      $read[9] = substr($read[9], 3);  
      $read[10] = substr($read[10], 3);  
  } elsif ($read[1] == 77 && $correctBCtrans ne $correctBC) {  
      ($ub_tag) = ($read[12] =~ /UB:Z:(.*?)$/);  
      ($qu_tag) = ($read[14] =~ /QU:Z:(.*?)$/); 
      if (defined $qu_tag) { 
        $read[9] = $ub_tag . $read[9];  
        $read[10] = $qu_tag . $read[10];  
        $read[12] = "UB:Z:";  
        $read[14] = "QU:Z:"; 
      }
  } elsif ($read[1] == 141 && $correctBCtrans ne $correctBC) {
      ($ub_tag) = ($read[12] =~ /UB:Z:(.*?)$/);  
      ($qu_tag) = ($read[14] =~ /QU:Z:(.*?)$/); 
      if (defined $qu_tag) {
        $read[12] = "UB:Z:";  
        $read[14] = "QU:Z:";
      }
  }  

  
  print BCBAM $read[0],"\t",$read[1],"\t",$read[2],"\t",$read[3],"\t",$read[4],"\t",$read[5],"\t",$read[6],"\t",$read[7],"\t",$read[8],"\t",
        $read[9],"\t",$read[10],"\t","BX:Z:",$thisBC,"\t","BC:Z:",$correctBCtrans,"\t",$read[12],"\t",$read[13],"\t",$read[14],"\n";

}
close INBAM;
close BCBAM;
