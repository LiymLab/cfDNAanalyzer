for file in sites_hg19/*; do
    filename=$(basename "$file")
    ~/projects/202402_cfDNAIntegratedTool/input/liftOver "$file" ~/projects/202402_cfDNAIntegratedTool/input/hg38ToHg19.over.chain.gz sites_hg19_txt/"$filename" sites_hg19_left/left_"$filename"
done