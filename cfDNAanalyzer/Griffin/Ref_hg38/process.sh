for file in sites/*; do
    filename=$(basename "$file")
    tail -n +2 "$file" | awk -F '\t' '{print $1 "\t" $2 "\t" $2+1}' > "sites_hg19/$filename"
done