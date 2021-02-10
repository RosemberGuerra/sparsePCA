# Download and processing data for DR application #
# Rosember Guerra-Urzola 
# created: 11-11-2020
# edited:

# Using this file, we can download the data set from the
# url: 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7329/matrix/GSE7329_series_matrix.txt.gz'
# name of data set: GSE7329_series_matrix.txt.gz
# The data is downloaded in as .gz format
# The data 'GSE7329_series_matrix.txt.gz' is save in the same directory as this file.

rm(list = ls(all.names = TRUE))

###  Download data set ###
Data_loacation= 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7329/matrix/GSE7329_series_matrix.txt.gz'
destfile= paste0(getwd(),"/GSE7329_series_matrix.txt.gz")
download.file(Data_loacation,destfile)

full_data=gzfile('GSE7329_series_matrix.txt.gz','rt')  
full_data=read.csv(full_data,header=F,sep = '\t')

### Column names from the data ###
col_names= c()
for (i in 736:751) {
  # col_names=c(col_names, as.character(full_data[i,1]),as.character(dat[i,2]))
  col_names=c(col_names,sapply(full_data[i,], as.character))
} 
View(col_names)

### Data neede for the application ###
matrix1= full_data[752:703647,]
matrix2= matrix(rep(0, prod(dim(matrix1))), ncol = dim(matrix1)[2])
for(i in 1:nrow(matrix1) ){
  matrix2[i,]= as.numeric(sapply(matrix1[i,],as.character))
}

Data_Autism=matrix(data = t(matrix2),byrow = TRUE, ncol = length(col_names))
colnames(Data_Autism)=col_names # adding names
IndexCol= Data_Autism[,1]
Data_Autism=Data_Autism[,-c(1,dim(Data_Autism)[2])] # removing first and last column 
Data_Autism=Data_Autism[,-c(12,13,27)] # Removing mislabel individuals 
dim(Data_Autism)

### Removing NA values ###
Data_Autism= t(Data_Autism)
Data_Autism=scale(Data_Autism, center = TRUE, scale = TRUE)

NA_Data_Autism= is.na(Data_Autism)
NA_index= colSums(NA_Data_Autism)!=0
sum(NA_index) # number of NA columns
which(NA_index!=0)
Data_Autism=Data_Autism[,-which(NA_index!=0)]
dim(Data_Autism)

Data_Autism=scale(Data_Autism, center = TRUE, scale = TRUE)
Data_Autism=as.data.frame(t(Data_Autism))

### saving data set ###
write.csv(Data_Autism, file = 'Data_Autism.csv',sep = ',')
