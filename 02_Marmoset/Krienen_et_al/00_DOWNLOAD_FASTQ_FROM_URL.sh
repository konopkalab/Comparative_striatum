# importing packages
library(httr)
library(XML)
library(dplyr)

url = "https://data.nemoarchive.org/biccn/grant/u01_feng/mccarroll/transcriptome/sncell/10Xv3/marmoset/raw/2019-04-26_BI005_marm027_Spencer_Put_rxn2/"

# making http request
resource = GET(url)

# parsing data to html format
parse = htmlParse(resource)

# scrapping all the href tags
links = xpathSApply(parse,path = "//a",xmlGetAttr,"href")
links_sub = links[grepl('Put', links)] 

# All animals
gsub('.*marm', 'marm', links_sub) %>% gsub('_.*', '', .) %>% unique

# Extract the sample
links_sub2 = links_sub[grepl('marm027', links_sub)]
links_sub2 = paste0(url, links_sub2)

write.table(links_sub2, "/home2/gkonop/project/00_FASTQ/MARMOSET_KRIENEN/PUTAMEN/TMP/links_putamen.txt", row.names = F, col.names = F, quote = F)


## IN BASH ##
cd /home2/gkonop/project/00_FASTQ/MARMOSET_KRIENEN/PUTAMEN/TMP

for i in `cat links_putamen.txt`
do
wget $i
done
