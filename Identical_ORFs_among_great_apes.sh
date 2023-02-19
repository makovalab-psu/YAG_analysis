#!/bin/bash
#SBATCH --job-name=marta
#SBATCH --output=marta-%j.out
#SBATCH --error=marta-%j.err
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=500G

cat All_BPY2_orf_nuc.fa  | python fasta_to_hash_4.py --canonical --report=sequence > BPY2_orf_nuc.signatures 

cat BPY2_orf_nuc.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > BPY2_orf_nuc_shared.signatures
                
  cat BPY2_orf_nuc_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
	> BPY2_orf_nuc_shared.fasta
      
cat BPY2_orf_nuc_shared.fasta | grep ">" | wc -l > count_BPY2_orf_nuc_shared.txt

cat All_CDY_orf_nuc.fa  | python fasta_to_hash_4.py --canonical --report=sequence > CDY_orf_nuc.signatures 

cat CDY_orf_nuc.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > CDY_orf_nuc_shared.signatures
                
  cat CDY_orf_nuc_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > CDY_orf_nuc_shared.fasta
      
cat CDY_orf_nuc_shared.fasta | grep ">" | wc -l > count_CDY_orf_nuc_shared.txt

cat All_DAZ_orf_nuc.fa  | python fasta_to_hash_4.py --canonical --report=sequence > DAZ_orf_nuc.signatures 

cat DAZ_orf_nuc.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > DAZ_orf_nuc_shared.signatures
                
  cat DAZ_orf_nuc_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > DAZ_orf_nuc_shared.fasta
      
cat DAZ_orf_nuc_shared.fasta | grep ">" | wc -l > count_DAZ_orf_nuc_shared.txt


cat All_HSFY_orf_nuc.fa  | python fasta_to_hash_4.py --canonical --report=sequence > HSFY_orf_nuc.signatures 

cat HSFY_orf_nuc.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > HSFY_orf_nuc_shared.signatures
                
  cat HSFY_orf_nuc_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > HSFY_orf_nuc_shared.fasta
      
cat HSFY_orf_nuc_shared.fasta | grep ">" | wc -l > count_HSFY_orf_nuc_shared.txt


cat All_PRY_orf_nuc.fa  | python fasta_to_hash_4.py --canonical --report=sequence > PRY_orf_nuc.signatures 

cat PRY_orf_nuc.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > PRY_orf_nuc_shared.signatures
                
  cat PRY_orf_nuc_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > PRY_orf_nuc_shared.fasta
      
cat PRY_orf_nuc_shared.fasta | grep ">" | wc -l > count_PRY_orf_nuc_shared.txt


cat All_RBMY_orf_nuc.fa  | python fasta_to_hash_4.py --canonical --report=sequence > RBMY_orf_nuc.signatures 

cat RBMY_orf_nuc.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > RBMY_orf_nuc_shared.signatures
                
  cat RBMY_orf_nuc_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > RBMY_orf_nuc_shared.fasta
      
cat RBMY_orf_nuc_shared.fasta | grep ">" | wc -l > count_RBMY_orf_nuc_shared.txt


cat All_TSPY_orf_nuc.fa  | python fasta_to_hash_4.py --canonical --report=sequence > TSPY_orf_nuc.signatures 

cat TSPY_orf_nuc.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > TSPY_orf_nuc_shared.signatures
                
  cat TSPY_orf_nuc_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > TSPY_orf_nuc_shared.fasta
      
cat TSPY_orf_nuc_shared.fasta | grep ">" | wc -l > count_TSPY_orf_nuc_shared.txt


cat All_VCY_orf_nuc.fa  | python fasta_to_hash_4.py --canonical --report=sequence > VCY_orf_nuc.signatures 

cat VCY_orf_nuc.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > VCY_orf_nuc_shared.signatures
                
  cat VCY_orf_nuc_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > VCY_orf_nuc_shared.fasta
      
cat VCY_orf_nuc_shared.fasta | grep ">" | wc -l > count_VCY_orf_nuc_shared.txt


cat All_XKRY_orf_nuc.fa  | python fasta_to_hash_4.py --canonical --report=sequence > XKRY_orf_nuc.signatures 

cat XKRY_orf_nuc.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > XKRY_orf_nuc_shared.signatures
                
  cat XKRY_orf_nuc_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > XKRY_orf_nuc_shared.fasta
      
cat XKRY_orf_nuc_shared.fasta | grep ">" | wc -l > count_XKRY_orf_nuc_shared.txt

cat All_BPY2_orf_prot.fa  | python fasta_to_hash_4.py --canonical --report=sequence > BPY2_orf_prot.signatures 

cat BPY2_orf_prot.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > BPY2_orf_prot_shared.signatures
                
  cat BPY2_orf_prot_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > BPY2_orf_prot_shared.fasta
      
cat BPY2_orf_prot_shared.fasta | grep ">" | wc -l > count_BPY2_orf_prot_shared.txt

cat All_CDY_orf_prot.fa  | python fasta_to_hash_4.py --canonical --report=sequence > CDY_orf_prot.signatures 

cat CDY_orf_prot.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > CDY_orf_prot_shared.signatures
                
  cat CDY_orf_prot_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > CDY_orf_prot_shared.fasta
      
cat CDY_orf_prot_shared.fasta | grep ">" | wc -l > count_CDY_orf_prot_shared.txt

cat All_DAZ_orf_prot.fa  | python fasta_to_hash_4.py --canonical --report=sequence > DAZ_orf_prot.signatures 

cat DAZ_orf_prot.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > DAZ_orf_prot_shared.signatures
                
  cat DAZ_orf_prot_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > DAZ_orf_prot_shared.fasta
      
cat DAZ_orf_prot_shared.fasta | grep ">" | wc -l > count_DAZ_orf_prot_shared.txt


cat All_HSFY_orf_prot.fa  | python fasta_to_hash_4.py --canonical --report=sequence > HSFY_orf_prot.signatures 

cat HSFY_orf_prot.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > HSFY_orf_prot_shared.signatures
                
  cat HSFY_orf_prot_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > HSFY_orf_prot_shared.fasta
      
cat HSFY_orf_prot_shared.fasta | grep ">" | wc -l > count_HSFY_orf_prot_shared.txt


cat All_PRY_orf_prot.fa  | python fasta_to_hash_4.py --canonical --report=sequence > PRY_orf_prot.signatures 

cat PRY_orf_prot.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > PRY_orf_prot_shared.signatures
                
  cat PRY_orf_prot_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > PRY_orf_prot_shared.fasta
      
cat PRY_orf_prot_shared.fasta | grep ">" | wc -l > count_PRY_orf_prot_shared.txt


cat All_RBMY_orf_prot.fa  | python fasta_to_hash_4.py --canonical --report=sequence > RBMY_orf_prot.signatures 

cat RBMY_orf_prot.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > RBMY_orf_prot_shared.signatures
                
  cat RBMY_orf_prot_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > RBMY_orf_prot_shared.fasta
      
cat RBMY_orf_prot_shared.fasta | grep ">" | wc -l > count_RBMY_orf_prot_shared.txt


cat All_TSPY_orf_prot.fa  | python fasta_to_hash_4.py --canonical --report=sequence > TSPY_orf_prot.signatures 

cat TSPY_orf_prot.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > TSPY_orf_prot_shared.signatures
                
  cat TSPY_orf_prot_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > TSPY_orf_prot_shared.fasta
      
cat TSPY_orf_prot_shared.fasta | grep ">" | wc -l > count_TSPY_orf_prot_shared.txt


cat All_VCY_orf_prot.fa  | python fasta_to_hash_4.py --canonical --report=sequence > VCY_orf_prot.signatures 

cat VCY_orf_prot.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > VCY_orf_prot_shared.signatures
                
  cat VCY_orf_prot_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > VCY_orf_prot_shared.fasta
      
cat VCY_orf_prot_shared.fasta | grep ">" | wc -l > count_VCY_orf_prot_shared.txt


cat All_XKRY_orf_prot.fa  | python fasta_to_hash_4.py --canonical --report=sequence > XKRY_orf_prot.signatures 

cat XKRY_orf_prot.signatures \
      | sort \
      | awk '{
                hashVal = $1;
                if (hashVal != prevHashVal)
                  {
                  if (count > 1) print matches;
                  count = 1;
                  matches = $0;
                  prevHashVal = hashVal;
                  }
                else
                  {
                  count++;
                  matches = matches"\n"$0;
                  }
                }
         END {
                if (count > 1) print matches;
                }' \
      > XKRY_orf_prot_shared.signatures
                
  cat XKRY_orf_prot_shared.signatures \
      | awk '{
                hashVal = $1;  header=$2; nts=$3;
                if (hashVal != prevHashVal) { print header; print nts; }
                prevHashVal = hashVal;
                }' \
  > XKRY_orf_prot_shared.fasta
      
cat XKRY_orf_prot_shared.fasta | grep ">" | wc -l > count_XKRY_orf_prot_shared.txt


