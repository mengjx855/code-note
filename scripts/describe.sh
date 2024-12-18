#!/usr/bin/bash

if [ $# -ne 2 ];then
    echo "$0 -c/r <profile>"
    echo -e  "  -c columns"
    echo -e  "  -r row"
fi

while getopts "crR" arg;
do
    case $arg in
        r)
            perl -e 'open I, "$ARGV[0]";$x = readline, I;chomp($x);@m=split(/\t/, $x);$n=@m-1;while(<>){chomp;$s=0;@l=split/\t/;foreach $k(@l[1..$#l]){$s+=$k};$ave=$s/$n;print "$l[0]\t$s\t$ave\t$n\n"}' $2
            ;;
        c)
            perl -e 'open I, "$ARGV[0]";$t=readline, I;chomp($t);@x=split(/\t/,$t);%h;while(<I>){chomp;@l=split/\t/;foreach $k(1..@l-1){$h{$x[$k]}+=$l[$k]}}; $line=$.-1; foreach $k(keys %h){$ave=$h{$k}/$line;print "$k\t$h{$k}\t$ave\t$line\n"}' $2
            ;;
        R)
            perl -e 'open I, "$ARGV[0]";$x = readline, I;chomp($x);@m=split(/\t/, $x);$n=@m-1;while(<>){chomp;$c=0;$s=0;@l=split/\t/;foreach $k(@l[1..$#l]){$s+=$k;$c++ if $k >0 };$ave=$s/$c;print "$l[0]\t$s\t$ave\t$c\n"}' $2
            ;;
        \?)
            echo "other"
            ;;
    esac
done

