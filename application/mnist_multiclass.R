# Testing multiclass probability estimation #
# Data : MNIST #
rm(list = ls())

library(keras)
##########
## Data ##
##########
# Load the Data 
mnist_lst <- dataset_mnist()

## Load train_x,y
train_x <- mnist_lst$train$x
train_y <- mnist_lst$train$y

## Load test_x,y
test_x <- mnist_lst$test$x
test_y <- mnist_lst$test$y

## Data Flattening (2D array --> 1D array)
train_x <- 3*array(as.numeric(train_x), dim = c(dim(train_x)[[1]], 784))
test_x <- 3*array(as.numeric(test_x), dim = c(dim(test_x)[[1]], 784))

# Sampling train set
set.seed(1)
index <- sample(dim(train_x)[1], size = 1000, replace=FALSE)
train.x <- train_x[index,]
train.y <- train_y[index]

# Sampling test set
set.seed(1)
index1 <- sample(dim(test_x)[1], size = 500, replace=FALSE)
test.x <- test_x[index1,]
test.y <- test_y[index1]


# Divide train set pairwise --------------------------------------------------------------------(TRAIN)
#### Extract index pairwise####
#### class 0 pair
trn.idx01 <- which(train.y==0 | train.y==1)
trn.idx02 <- which(train.y==0 | train.y==2)
trn.idx03 <- which(train.y==0 | train.y==3)
trn.idx04 <- which(train.y==0 | train.y==4)
trn.idx05 <- which(train.y==0 | train.y==5)
trn.idx06 <- which(train.y==0 | train.y==6)
trn.idx07 <- which(train.y==0 | train.y==7)
trn.idx08 <- which(train.y==0 | train.y==8)
trn.idx09 <- which(train.y==0 | train.y==9)

#### class 1 pair
trn.idx12 <- which(train.y==1 | train.y==2)
trn.idx13 <- which(train.y==1 | train.y==3)
trn.idx14 <- which(train.y==1 | train.y==4)
trn.idx15 <- which(train.y==1 | train.y==5)
trn.idx16 <- which(train.y==1 | train.y==6)
trn.idx17 <- which(train.y==1 | train.y==7)
trn.idx18 <- which(train.y==1 | train.y==8)
trn.idx19 <- which(train.y==1 | train.y==9)

#### class 2 pair
trn.idx23 <- which(train.y==2 | train.y==3)
trn.idx24 <- which(train.y==2 | train.y==4)
trn.idx25 <- which(train.y==2 | train.y==5)
trn.idx26 <- which(train.y==2 | train.y==6)
trn.idx27 <- which(train.y==2 | train.y==7)
trn.idx28 <- which(train.y==2 | train.y==8)
trn.idx29 <- which(train.y==2 | train.y==9)

#### class 3 pair
trn.idx34 <- which(train.y==3 | train.y==4)
trn.idx35 <- which(train.y==3 | train.y==5)
trn.idx36 <- which(train.y==3 | train.y==6)
trn.idx37 <- which(train.y==3 | train.y==7)
trn.idx38 <- which(train.y==3 | train.y==8)
trn.idx39 <- which(train.y==3 | train.y==9)

#### class 4 pair
trn.idx45 <- which(train.y==4 | train.y==5)
trn.idx46 <- which(train.y==4 | train.y==6)
trn.idx47 <- which(train.y==4 | train.y==7)
trn.idx48 <- which(train.y==4 | train.y==8)
trn.idx49 <- which(train.y==4 | train.y==9)

#### class 5 pair
trn.idx56 <- which(train.y==5 | train.y==6)
trn.idx57 <- which(train.y==5 | train.y==7)
trn.idx58 <- which(train.y==5 | train.y==8)
trn.idx59 <- which(train.y==5 | train.y==9)

#### class 6 pair
trn.idx67 <- which(train.y==6 | train.y==7)
trn.idx68 <- which(train.y==6 | train.y==8)
trn.idx69 <- which(train.y==6 | train.y==9)

#### class 7 pair
trn.idx78 <- which(train.y==7 | train.y==8)
trn.idx79 <- which(train.y==7 | train.y==9)

#### class 8 pair
trn.idx89 <- which(train.y==8 | train.y==9)

#### Split train.x pariwise ####
#### class 0 pair
train.x01 <- train.x[trn.idx01,]
train.x02 <- train.x[trn.idx02,]
train.x03 <- train.x[trn.idx03,]
train.x04 <- train.x[trn.idx04,]
train.x05 <- train.x[trn.idx05,]
train.x06 <- train.x[trn.idx06,]
train.x07 <- train.x[trn.idx07,]
train.x08 <- train.x[trn.idx08,]
train.x09 <- train.x[trn.idx09,]
#### class 1 pair
train.x12 <- train.x[trn.idx12,]
train.x13 <- train.x[trn.idx13,]
train.x14 <- train.x[trn.idx14,]
train.x15 <- train.x[trn.idx15,]
train.x16 <- train.x[trn.idx16,]
train.x17 <- train.x[trn.idx17,]
train.x18 <- train.x[trn.idx18,]
train.x19 <- train.x[trn.idx19,]
#### class 2 pair
train.x23 <- train.x[trn.idx23,]
train.x24 <- train.x[trn.idx24,]
train.x25 <- train.x[trn.idx25,]
train.x26 <- train.x[trn.idx26,]
train.x27 <- train.x[trn.idx27,]
train.x28 <- train.x[trn.idx28,]
train.x29 <- train.x[trn.idx29,]
#### class 3 pair
train.x34 <- train.x[trn.idx34,]
train.x35 <- train.x[trn.idx35,]
train.x36 <- train.x[trn.idx36,]
train.x37 <- train.x[trn.idx37,]
train.x38 <- train.x[trn.idx38,]
train.x39 <- train.x[trn.idx39,]
#### class 4 pair
train.x45 <- train.x[trn.idx45,]
train.x46 <- train.x[trn.idx46,]
train.x47 <- train.x[trn.idx47,]
train.x48 <- train.x[trn.idx48,]
train.x49 <- train.x[trn.idx49,]
#### class 5 pair
train.x56 <- train.x[trn.idx56,]
train.x57 <- train.x[trn.idx57,]
train.x58 <- train.x[trn.idx58,]
train.x59 <- train.x[trn.idx59,]
#### class 6 pair
train.x67 <- train.x[trn.idx67,]
train.x68 <- train.x[trn.idx68,]
train.x69 <- train.x[trn.idx69,]
#### class 7 pair
train.x78 <- train.x[trn.idx78,]
train.x79 <- train.x[trn.idx79,]
#### class 8 pair
train.x89 <- train.x[trn.idx89,]

#### Split train.y pariwise ####
#### class 0 pair
train.y01 <- train.y[trn.idx01]
train.y02 <- train.y[trn.idx02]
train.y03 <- train.y[trn.idx03]
train.y04 <- train.y[trn.idx04]
train.y05 <- train.y[trn.idx05]
train.y06 <- train.y[trn.idx06]
train.y07 <- train.y[trn.idx07]
train.y08 <- train.y[trn.idx08]
train.y09 <- train.y[trn.idx09]
#### class 1 pair
train.y12 <- train.y[trn.idx12]
train.y13 <- train.y[trn.idx13]
train.y14 <- train.y[trn.idx14]
train.y15 <- train.y[trn.idx15]
train.y16 <- train.y[trn.idx16]
train.y17 <- train.y[trn.idx17]
train.y18 <- train.y[trn.idx18]
train.y19 <- train.y[trn.idx19]
#### class 2 pair
train.y23 <- train.y[trn.idx23]
train.y24 <- train.y[trn.idx24]
train.y25 <- train.y[trn.idx25]
train.y26 <- train.y[trn.idx26]
train.y27 <- train.y[trn.idx27]
train.y28 <- train.y[trn.idx28]
train.y29 <- train.y[trn.idx29]
#### class 3 pair
train.y34 <- train.y[trn.idx34]
train.y35 <- train.y[trn.idx35]
train.y36 <- train.y[trn.idx36]
train.y37 <- train.y[trn.idx37]
train.y38 <- train.y[trn.idx38]
train.y39 <- train.y[trn.idx39]
#### class 4 pair
train.y45 <- train.y[trn.idx45]
train.y46 <- train.y[trn.idx46]
train.y47 <- train.y[trn.idx47]
train.y48 <- train.y[trn.idx48]
train.y49 <- train.y[trn.idx49]
#### class 5 pair
train.y56 <- train.y[trn.idx56]
train.y57 <- train.y[trn.idx57]
train.y58 <- train.y[trn.idx58]
train.y59 <- train.y[trn.idx59]
#### class 6 pair
train.y67 <- train.y[trn.idx67]
train.y68 <- train.y[trn.idx68]
train.y69 <- train.y[trn.idx69]
#### class 7 pair
train.y78 <- train.y[trn.idx78]
train.y79 <- train.y[trn.idx79]
#### class 8 pair
train.y89 <- train.y[trn.idx89]

#### Change y value into -1 and 1 ####
#### clas 0 pair
train.y01 <- ifelse(train.y01 == 1, 1, -1)
train.y02 <- ifelse(train.y02 == 2, 1, -1)
train.y03 <- ifelse(train.y03 == 3, 1, -1)
train.y04 <- ifelse(train.y04 == 4, 1, -1)
train.y05 <- ifelse(train.y05 == 5, 1, -1)
train.y06 <- ifelse(train.y06 == 6, 1, -1)
train.y07 <- ifelse(train.y07 == 7, 1, -1)
train.y08 <- ifelse(train.y08 == 8, 1, -1)
train.y09 <- ifelse(train.y09 == 9, 1, -1)
#### clas 1 pair
train.y12 <- ifelse(train.y12 == 2, 1, -1)
train.y13 <- ifelse(train.y13 == 3, 1, -1)
train.y14 <- ifelse(train.y14 == 4, 1, -1)
train.y15 <- ifelse(train.y15 == 5, 1, -1)
train.y16 <- ifelse(train.y16 == 6, 1, -1)
train.y17 <- ifelse(train.y17 == 7, 1, -1)
train.y18 <- ifelse(train.y18 == 8, 1, -1)
train.y19 <- ifelse(train.y19 == 9, 1, -1)
#### clas 2 pair
train.y23 <- ifelse(train.y23 == 3, 1, -1)
train.y24 <- ifelse(train.y24 == 4, 1, -1)
train.y25 <- ifelse(train.y25 == 5, 1, -1)
train.y26 <- ifelse(train.y26 == 6, 1, -1)
train.y27 <- ifelse(train.y27 == 7, 1, -1)
train.y28 <- ifelse(train.y28 == 8, 1, -1)
train.y29 <- ifelse(train.y29 == 9, 1, -1)
#### clas 3 pair
train.y34 <- ifelse(train.y34 == 4, 1, -1)
train.y35 <- ifelse(train.y35 == 5, 1, -1)
train.y36 <- ifelse(train.y36 == 6, 1, -1)
train.y37 <- ifelse(train.y37 == 7, 1, -1)
train.y38 <- ifelse(train.y38 == 8, 1, -1)
train.y39 <- ifelse(train.y39 == 9, 1, -1)
#### clas 4 pair
train.y45 <- ifelse(train.y45 == 5, 1, -1)
train.y46 <- ifelse(train.y46 == 6, 1, -1)
train.y47 <- ifelse(train.y47 == 7, 1, -1)
train.y48 <- ifelse(train.y48 == 8, 1, -1)
train.y49 <- ifelse(train.y49 == 9, 1, -1)
#### clas 5 pair
train.y56 <- ifelse(train.y56 == 6, 1, -1)
train.y57 <- ifelse(train.y57 == 7, 1, -1)
train.y58 <- ifelse(train.y58 == 8, 1, -1)
train.y59 <- ifelse(train.y59 == 9, 1, -1)
#### clas 6 pair
train.y67 <- ifelse(train.y67 == 7, 1, -1)
train.y68 <- ifelse(train.y68 == 8, 1, -1)
train.y69 <- ifelse(train.y69 == 9, 1, -1)
#### clas 7 pair
train.y78 <- ifelse(train.y78 == 8, 1, -1)
train.y79 <- ifelse(train.y79 == 9, 1, -1)
#### clas 8 pair
train.y89 <- ifelse(train.y89 == 9, 1, -1)



# Divide TEST set pairwise ----------------------------------------------------------------------(TEST)
#### Extract index pairwise####
#### class 0 pair
tst.idx01 <- which(test.y==0 | test.y==1)
tst.idx02 <- which(test.y==0 | test.y==2)
tst.idx03 <- which(test.y==0 | test.y==3)
tst.idx04 <- which(test.y==0 | test.y==4)
tst.idx05 <- which(test.y==0 | test.y==5)
tst.idx06 <- which(test.y==0 | test.y==6)
tst.idx07 <- which(test.y==0 | test.y==7)
tst.idx08 <- which(test.y==0 | test.y==8)
tst.idx09 <- which(test.y==0 | test.y==9)
#### class 1 pair
tst.idx12 <- which(test.y==1 | test.y==2)
tst.idx13 <- which(test.y==1 | test.y==3)
tst.idx14 <- which(test.y==1 | test.y==4)
tst.idx15 <- which(test.y==1 | test.y==5)
tst.idx16 <- which(test.y==1 | test.y==6)
tst.idx17 <- which(test.y==1 | test.y==7)
tst.idx18 <- which(test.y==1 | test.y==8)
tst.idx19 <- which(test.y==1 | test.y==9)
#### class 2 pair
tst.idx23 <- which(test.y==2 | test.y==3)
tst.idx24 <- which(test.y==2 | test.y==4)
tst.idx25 <- which(test.y==2 | test.y==5)
tst.idx26 <- which(test.y==2 | test.y==6)
tst.idx27 <- which(test.y==2 | test.y==7)
tst.idx28 <- which(test.y==2 | test.y==8)
tst.idx29 <- which(test.y==2 | test.y==9)
#### class 3 pair
tst.idx34 <- which(test.y==3 | test.y==4)
tst.idx35 <- which(test.y==3 | test.y==5)
tst.idx36 <- which(test.y==3 | test.y==6)
tst.idx37 <- which(test.y==3 | test.y==7)
tst.idx38 <- which(test.y==3 | test.y==8)
tst.idx39 <- which(test.y==3 | test.y==9)
#### class 4 pair
tst.idx45 <- which(test.y==4 | test.y==5)
tst.idx46 <- which(test.y==4 | test.y==6)
tst.idx47 <- which(test.y==4 | test.y==7)
tst.idx48 <- which(test.y==4 | test.y==8)
tst.idx49 <- which(test.y==4 | test.y==9)
#### class 5 pair
tst.idx56 <- which(test.y==5 | test.y==6)
tst.idx57 <- which(test.y==5 | test.y==7)
tst.idx58 <- which(test.y==5 | test.y==8)
tst.idx59 <- which(test.y==5 | test.y==9)
#### class 6 pair
tst.idx67 <- which(test.y==6 | test.y==7)
tst.idx68 <- which(test.y==6 | test.y==8)
tst.idx69 <- which(test.y==6 | test.y==9)
#### class 7 pair
tst.idx78 <- which(test.y==7 | test.y==8)
tst.idx79 <- which(test.y==7 | test.y==9)
#### class 8 pair
tst.idx89 <- which(test.y==8 | test.y==9)

#### Split test.x pariwise ####
#### class 0 pair
test.x01 <- test.x[tst.idx01]
test.x02 <- test.x[tst.idx02]
test.x03 <- test.x[tst.idx03]
test.x04 <- test.x[tst.idx04]
test.x05 <- test.x[tst.idx05]
test.x06 <- test.x[tst.idx06]
test.x07 <- test.x[tst.idx07]
test.x08 <- test.x[tst.idx08]
test.x09 <- test.x[tst.idx09]
#### class 1 pair
test.x12 <- test.x[tst.idx12]
test.x13 <- test.x[tst.idx13]
test.x14 <- test.x[tst.idx14]
test.x15 <- test.x[tst.idx15]
test.x16 <- test.x[tst.idx16]
test.x17 <- test.x[tst.idx17]
test.x18 <- test.x[tst.idx18]
test.x19 <- test.x[tst.idx19]
#### class 2 pair
test.x23 <- test.x[tst.idx23]
test.x24 <- test.x[tst.idx24]
test.x25 <- test.x[tst.idx25]
test.x26 <- test.x[tst.idx26]
test.x27 <- test.x[tst.idx27]
test.x28 <- test.x[tst.idx28]
test.x29 <- test.x[tst.idx29]
#### class 3 pair
test.x34 <- test.x[tst.idx34]
test.x35 <- test.x[tst.idx35]
test.x36 <- test.x[tst.idx36]
test.x37 <- test.x[tst.idx37]
test.x38 <- test.x[tst.idx38]
test.x39 <- test.x[tst.idx39]
#### class 4 pair
test.x45 <- test.x[tst.idx45]
test.x46 <- test.x[tst.idx46]
test.x47 <- test.x[tst.idx47]
test.x48 <- test.x[tst.idx48]
test.x49 <- test.x[tst.idx49]
#### class 5 pair
test.x56 <- test.x[tst.idx56]
test.x57 <- test.x[tst.idx57]
test.x58 <- test.x[tst.idx58]
test.x59 <- test.x[tst.idx59]
#### class 6 pair
test.x67 <- test.x[tst.idx67]
test.x68 <- test.x[tst.idx68]
test.x69 <- test.x[tst.idx69]
#### class 7 pair
test.x78 <- test.x[tst.idx78]
test.x79 <- test.x[tst.idx79]
#### class 8 pair
test.x89 <- test.x[tst.idx89]


#### Split test.y pariwise ####
#### class 0 pair
test.y01 <- test.y[tst.idx01]
test.y02 <- test.y[tst.idx02]
test.y03 <- test.y[tst.idx03]
test.y04 <- test.y[tst.idx04]
test.y05 <- test.y[tst.idx05]
test.y06 <- test.y[tst.idx06]
test.y07 <- test.y[tst.idx07]
test.y08 <- test.y[tst.idx08]
test.y09 <- test.y[tst.idx09]
#### class 1 pair
test.y12 <- test.y[tst.idx12]
test.y13 <- test.y[tst.idx13]
test.y14 <- test.y[tst.idx14]
test.y15 <- test.y[tst.idx15]
test.y16 <- test.y[tst.idx16]
test.y17 <- test.y[tst.idx17]
test.y18 <- test.y[tst.idx18]
test.y19 <- test.y[tst.idx19]
#### class 2 pair
test.y23 <- test.y[tst.idx23]
test.y24 <- test.y[tst.idx24]
test.y25 <- test.y[tst.idx25]
test.y26 <- test.y[tst.idx26]
test.y27 <- test.y[tst.idx27]
test.y28 <- test.y[tst.idx28]
test.y29 <- test.y[tst.idx29]
#### class 3 pair
test.y34 <- test.y[tst.idx34]
test.y35 <- test.y[tst.idx35]
test.y36 <- test.y[tst.idx36]
test.y37 <- test.y[tst.idx37]
test.y38 <- test.y[tst.idx38]
test.y39 <- test.y[tst.idx39]
#### class 4 pair
test.y45 <- test.y[tst.idx45]
test.y46 <- test.y[tst.idx46]
test.y47 <- test.y[tst.idx47]
test.y48 <- test.y[tst.idx48]
test.y49 <- test.y[tst.idx49]
#### class 5 pair
test.y56 <- test.y[tst.idx56]
test.y57 <- test.y[tst.idx57]
test.y58 <- test.y[tst.idx58]
test.y59 <- test.y[tst.idx59]
#### class 6 pair
test.y67 <- test.y[tst.idx67]
test.y68 <- test.y[tst.idx68]
test.y69 <- test.y[tst.idx69]
#### class 7 pair
test.y78 <- test.y[tst.idx78]
test.y79 <- test.y[tst.idx79]
#### class 8 pair
test.y89 <- test.y[tst.idx89]

#### Change y value into -1 and 1 ####
#### clas 0 pair
test.y01 <- ifelse(test.y01 == 1, 1, -1)
test.y02 <- ifelse(test.y02 == 2, 1, -1)
test.y03 <- ifelse(test.y03 == 3, 1, -1)
test.y04 <- ifelse(test.y01 == 4, 1, -1)
test.y05 <- ifelse(test.y02 == 5, 1, -1)
test.y06 <- ifelse(test.y03 == 6, 1, -1)
test.y07 <- ifelse(test.y01 == 7, 1, -1)
test.y08 <- ifelse(test.y02 == 8, 1, -1)
test.y09 <- ifelse(test.y03 == 9, 1, -1)
#### clas 1 pair
test.y12 <- ifelse(test.y12 == 2, 1, -1)
test.y13 <- ifelse(test.y13 == 3, 1, -1)
test.y14 <- ifelse(test.y14 == 4, 1, -1)
test.y15 <- ifelse(test.y15 == 5, 1, -1)
test.y16 <- ifelse(test.y16 == 6, 1, -1)
test.y17 <- ifelse(test.y17 == 7, 1, -1)
test.y18 <- ifelse(test.y18 == 8, 1, -1)
test.y19 <- ifelse(test.y19 == 9, 1, -1)
#### clas 2 pair
test.y23 <- ifelse(test.y23 == 3, 1, -1)
test.y24 <- ifelse(test.y24 == 4, 1, -1)
test.y25 <- ifelse(test.y25 == 5, 1, -1)
test.y26 <- ifelse(test.y26 == 6, 1, -1)
test.y27 <- ifelse(test.y27 == 7, 1, -1)
test.y28 <- ifelse(test.y28 == 8, 1, -1)
test.y29 <- ifelse(test.y29 == 9, 1, -1)
#### clas 3 pair
test.y34 <- ifelse(test.y34 == 4, 1, -1)
test.y35 <- ifelse(test.y35 == 5, 1, -1)
test.y36 <- ifelse(test.y36 == 6, 1, -1)
test.y37 <- ifelse(test.y37 == 7, 1, -1)
test.y38 <- ifelse(test.y38 == 8, 1, -1)
test.y39 <- ifelse(test.y39 == 9, 1, -1)
#### clas 4 pair
test.y45 <- ifelse(test.y45 == 5, 1, -1)
test.y46 <- ifelse(test.y46 == 6, 1, -1)
test.y47 <- ifelse(test.y47 == 7, 1, -1)
test.y48 <- ifelse(test.y48 == 8, 1, -1)
test.y49 <- ifelse(test.y49 == 9, 1, -1)
#### clas 5 pair
test.y56 <- ifelse(test.y56 == 6, 1, -1)
test.y57 <- ifelse(test.y57 == 7, 1, -1)
test.y58 <- ifelse(test.y58 == 8, 1, -1)
test.y59 <- ifelse(test.y59 == 9, 1, -1)
#### clas 6 pair
test.y67 <- ifelse(test.y67 == 7, 1, -1)
test.y68 <- ifelse(test.y68 == 8, 1, -1)
test.y69 <- ifelse(test.y69 == 9, 1, -1)
#### clas 7 pair
test.y78 <- ifelse(test.y78 == 8, 1, -1)
test.y79 <- ifelse(test.y79 == 9, 1, -1)
#### clas 8 pair
test.y89 <- ifelse(test.y89 == 9, 1, -1)


######################
## load the methods ##
######################
library(fda)
setwd('C:/Users/eugene/Desktop/SVM/shared/R code/')
source('eu/fsvm.prob.R')
source('eu/predict.fsvm.prob.R')
source('fn/fsvm.pi.path.R')
source('fn/fsvm.sub.pi.path.R')
dyn.load("C:/Users/eugene/Desktop/SVM_R/shared/R code/KernSurf/temp/wsvmqp.dll")
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('C:/Users/eugene/Desktop/SVM_R/shared/R code/KernSurf/R')

# Parameters
L <- 10
t <- seq(0, 1, by = 1/783)

####=======================     train data      ===============================####
#### convert matrix into train.x.list ####
# class 0 pair
num <- dim(train.x01)[1]
train.x01.list <- as.list(1:num)
for (i in 1:num) train.x01.list[[i]] <- matrix(train.x01[i,], nrow=1,byrow=TRUE)

num <- dim(train.x02)[1]
train.x02.list <- as.list(1:num)
for (i in 1:num) train.x02.list[[i]] <- matrix(train.x02[i,], nrow=1,byrow=TRUE)

num <- dim(train.x03)[1]
train.x03.list <- as.list(1:num)
for (i in 1:num) train.x03.list[[i]] <- matrix(train.x03[i,], nrow=1,byrow=TRUE)

num <- dim(train.x04)[1]
train.x04.list <- as.list(1:num)
for (i in 1:num) train.x04.list[[i]] <- matrix(train.x04[i,], nrow=1,byrow=TRUE)

num <- dim(train.x05)[1]
train.x05.list <- as.list(1:num)
for (i in 1:num) train.x05.list[[i]] <- matrix(train.x05[i,], nrow=1,byrow=TRUE)

num <- dim(train.x06)[1]
train.x06.list <- as.list(1:num)
for (i in 1:num) train.x06.list[[i]] <- matrix(train.x06[i,], nrow=1,byrow=TRUE)

num <- dim(train.x07)[1]
train.x07.list <- as.list(1:num)
for (i in 1:num) train.x07.list[[i]] <- matrix(train.x07[i,], nrow=1,byrow=TRUE)

num <- dim(train.x08)[1]
train.x08.list <- as.list(1:num)
for (i in 1:num) train.x08.list[[i]] <- matrix(train.x08[i,], nrow=1,byrow=TRUE)

num <- dim(train.x09)[1]
train.x09.list <- as.list(1:num)
for (i in 1:num) train.x09.list[[i]] <- matrix(train.x09[i,], nrow=1,byrow=TRUE)

# class 1 pair
num <- dim(train.x12)[1]
train.x12.list <- as.list(1:num)
for (i in 1:num) train.x12.list[[i]] <- matrix(train.x12[i,], nrow=1,byrow=TRUE)

num <- dim(train.x13)[1]
train.x13.list <- as.list(1:num)
for (i in 1:num) train.x13.list[[i]] <- matrix(train.x13[i,], nrow=1,byrow=TRUE)

num <- dim(train.x14)[1]
train.x14.list <- as.list(1:num)
for (i in 1:num) train.x14.list[[i]] <- matrix(train.x14[i,], nrow=1,byrow=TRUE)

num <- dim(train.x15)[1]
train.x15.list <- as.list(1:num)
for (i in 1:num) train.x15.list[[i]] <- matrix(train.x15[i,], nrow=1,byrow=TRUE)

num <- dim(train.x16)[1]
train.x16.list <- as.list(1:num)
for (i in 1:num) train.x16.list[[i]] <- matrix(train.x16[i,], nrow=1,byrow=TRUE)

num <- dim(train.x17)[1]
train.x17.list <- as.list(1:num)
for (i in 1:num) train.x17.list[[i]] <- matrix(train.x17[i,], nrow=1,byrow=TRUE)

num <- dim(train.x18)[1]
train.x18.list <- as.list(1:num)
for (i in 1:num) train.x18.list[[i]] <- matrix(train.x18[i,], nrow=1,byrow=TRUE)

num <- dim(train.x19)[1]
train.x19.list <- as.list(1:num)
for (i in 1:num) train.x19.list[[i]] <- matrix(train.x19[i,], nrow=1,byrow=TRUE)

# class 2 pair
num <- dim(train.x23)[1]
train.x23.list <- as.list(1:num)
for (i in 1:num) train.x23.list[[i]] <- matrix(train.x23[i,], nrow=1,byrow=TRUE)
num <- dim(train.x24)[1]
train.x24.list <- as.list(1:num)
for (i in 1:num) train.x24.list[[i]] <- matrix(train.x24[i,], nrow=1,byrow=TRUE)
num <- dim(train.x25)[1]
train.x25.list <- as.list(1:num)
for (i in 1:num) train.x25.list[[i]] <- matrix(train.x25[i,], nrow=1,byrow=TRUE)
num <- dim(train.x26)[1]
train.x26.list <- as.list(1:num)
for (i in 1:num) train.x26.list[[i]] <- matrix(train.x26[i,], nrow=1,byrow=TRUE)
num <- dim(train.x27)[1]
train.x27.list <- as.list(1:num)
for (i in 1:num) train.x27.list[[i]] <- matrix(train.x27[i,], nrow=1,byrow=TRUE)
num <- dim(train.x28)[1]
train.x28.list <- as.list(1:num)
for (i in 1:num) train.x28.list[[i]] <- matrix(train.x28[i,], nrow=1,byrow=TRUE)
num <- dim(train.x29)[1]
train.x29.list <- as.list(1:num)
for (i in 1:num) train.x29.list[[i]] <- matrix(train.x29[i,], nrow=1,byrow=TRUE)
# class 3 pair
num <- dim(train.x34)[1]
train.x34.list <- as.list(1:num)
for (i in 1:num) train.x34.list[[i]] <- matrix(train.x34[i,], nrow=1,byrow=TRUE)
num <- dim(train.x35)[1]
train.x35.list <- as.list(1:num)
for (i in 1:num) train.x35.list[[i]] <- matrix(train.x35[i,], nrow=1,byrow=TRUE)
num <- dim(train.x36)[1]
train.x36.list <- as.list(1:num)
for (i in 1:num) train.x36.list[[i]] <- matrix(train.x36[i,], nrow=1,byrow=TRUE)
num <- dim(train.x37)[1]
train.x37.list <- as.list(1:num)
for (i in 1:num) train.x37.list[[i]] <- matrix(train.x37[i,], nrow=1,byrow=TRUE)
num <- dim(train.x38)[1]
train.x38.list <- as.list(1:num)
for (i in 1:num) train.x38.list[[i]] <- matrix(train.x38[i,], nrow=1,byrow=TRUE)
num <- dim(train.x39)[1]
train.x39.list <- as.list(1:num)
for (i in 1:num) train.x39.list[[i]] <- matrix(train.x39[i,], nrow=1,byrow=TRUE)
# class 4 pair
num <- dim(train.x45)[1]
train.x45.list <- as.list(1:num)
for (i in 1:num) train.x45.list[[i]] <- matrix(train.x45[i,], nrow=1,byrow=TRUE)
num <- dim(train.x46)[1]
train.x46.list <- as.list(1:num)
for (i in 1:num) train.x46.list[[i]] <- matrix(train.x46[i,], nrow=1,byrow=TRUE)
num <- dim(train.x47)[1]
train.x47.list <- as.list(1:num)
for (i in 1:num) train.x47.list[[i]] <- matrix(train.x47[i,], nrow=1,byrow=TRUE)
num <- dim(train.x48)[1]
train.x48.list <- as.list(1:num)
for (i in 1:num) train.x48.list[[i]] <- matrix(train.x48[i,], nrow=1,byrow=TRUE)
num <- dim(train.x49)[1]
train.x49.list <- as.list(1:num)
for (i in 1:num) train.x49.list[[i]] <- matrix(train.x49[i,], nrow=1,byrow=TRUE)
# class 5 pair
num <- dim(train.x56)[1]
train.x56.list <- as.list(1:num)
for (i in 1:num) train.x56.list[[i]] <- matrix(train.x56[i,], nrow=1,byrow=TRUE)
num <- dim(train.x57)[1]
train.x57.list <- as.list(1:num)
for (i in 1:num) train.x57.list[[i]] <- matrix(train.x57[i,], nrow=1,byrow=TRUE)
num <- dim(train.x58)[1]
train.x58.list <- as.list(1:num)
for (i in 1:num) train.x58.list[[i]] <- matrix(train.x58[i,], nrow=1,byrow=TRUE)
num <- dim(train.x59)[1]
train.x59.list <- as.list(1:num)
for (i in 1:num) train.x59.list[[i]] <- matrix(train.x59[i,], nrow=1,byrow=TRUE)
# class 6 pair
num <- dim(train.x67)[1]
train.x67.list <- as.list(1:num)
for (i in 1:num) train.x67.list[[i]] <- matrix(train.x67[i,], nrow=1,byrow=TRUE)
num <- dim(train.x68)[1]
train.x68.list <- as.list(1:num)
for (i in 1:num) train.x68.list[[i]] <- matrix(train.x68[i,], nrow=1,byrow=TRUE)
num <- dim(train.x69)[1]
train.x69.list <- as.list(1:num)
for (i in 1:num) train.x69.list[[i]] <- matrix(train.x69[i,], nrow=1,byrow=TRUE)
# class 7 pair
num <- dim(train.x78)[1]
train.x78.list <- as.list(1:num)
for (i in 1:num) train.x78.list[[i]] <- matrix(train.x78[i,], nrow=1,byrow=TRUE)
num <- dim(train.x79)[1]
train.x79.list <- as.list(1:num)
for (i in 1:num) train.x79.list[[i]] <- matrix(train.x79[i,], nrow=1,byrow=TRUE)
# class 8 pair
num <- dim(train.x89)[1]
train.x89.list <- as.list(1:num)
for (i in 1:num) train.x89.list[[i]] <- matrix(train.x89[i,], nrow=1,byrow=TRUE)

#### calculate the pi path ####
# class 0 pair
obj01 <- fsvm.prob(train.x01.list, c(train.y01), t, L)
obj02 <- fsvm.prob(train.x02.list, c(train.y02), t, L)
obj03 <- fsvm.prob(train.x03.list, c(train.y03), t, L)
obj04 <- fsvm.prob(train.x04.list, c(train.y04), t, L)
obj05 <- fsvm.prob(train.x05.list, c(train.y05), t, L)
obj06 <- fsvm.prob(train.x06.list, c(train.y06), t, L)
obj07 <- fsvm.prob(train.x07.list, c(train.y07), t, L)
obj08 <- fsvm.prob(train.x08.list, c(train.y08), t, L)
obj09 <- fsvm.prob(train.x09.list, c(train.y09), t, L)
# class 1 pair
obj12 <- fsvm.prob(train.x12.list, c(train.y12), t, L)
obj13 <- fsvm.prob(train.x13.list, c(train.y13), t, L) #########################
obj14 <- fsvm.prob(train.x14.list, c(train.y14), t, L)
obj15 <- fsvm.prob(train.x15.list, c(train.y15), t, L)
obj16 <- fsvm.prob(train.x16.list, c(train.y16), t, L)
obj17 <- fsvm.prob(train.x17.list, c(train.y17), t, L)
obj18 <- fsvm.prob(train.x18.list, c(train.y18), t, L)
obj19 <- fsvm.prob(train.x19.list, c(train.y19), t, L)
# class 2 pair
obj23 <- fsvm.prob(train.x23.list, c(train.y23), t, L)
obj24 <- fsvm.prob(train.x24.list, c(train.y24), t, L)
obj25 <- fsvm.prob(train.x25.list, c(train.y25), t, L)
obj26 <- fsvm.prob(train.x26.list, c(train.y26), t, L)
obj27 <- fsvm.prob(train.x27.list, c(train.y27), t, L)
obj28 <- fsvm.prob(train.x28.list, c(train.y28), t, L)
obj29 <- fsvm.prob(train.x29.list, c(train.y29), t, L)
# class 3 pair
obj34 <- fsvm.prob(train.x34.list, c(train.y34), t, L)
obj35 <- fsvm.prob(train.x35.list, c(train.y35), t, L)
obj36 <- fsvm.prob(train.x36.list, c(train.y36), t, L)
obj37 <- fsvm.prob(train.x37.list, c(train.y37), t, L)
obj38 <- fsvm.prob(train.x38.list, c(train.y38), t, L)
obj39 <- fsvm.prob(train.x39.list, c(train.y39), t, L)
# class 4 pair
obj45 <- fsvm.prob(train.x45.list, c(train.y45), t, L)
obj46 <- fsvm.prob(train.x46.list, c(train.y46), t, L)
obj47 <- fsvm.prob(train.x47.list, c(train.y47), t, L)
obj48 <- fsvm.prob(train.x48.list, c(train.y48), t, L)
obj49 <- fsvm.prob(train.x49.list, c(train.y49), t, L)
# class 5 pair
obj56 <- fsvm.prob(train.x56.list, c(train.y56), t, L)
obj57 <- fsvm.prob(train.x57.list, c(train.y57), t, L)
obj58 <- fsvm.prob(train.x58.list, c(train.y58), t, L)
obj59 <- fsvm.prob(train.x59.list, c(train.y59), t, L)
# class 6 pair
obj67 <- fsvm.prob(train.x67.list, c(train.y67), t, L)
obj68 <- fsvm.prob(train.x68.list, c(train.y68), t, L)
obj69 <- fsvm.prob(train.x69.list, c(train.y69), t, L)
# class 7 pair
obj78 <- fsvm.prob(train.x78.list, c(train.y78), t, L)
obj79 <- fsvm.prob(train.x79.list, c(train.y79), t, L)
# class 8 pair
obj89 <- fsvm.prob(train.x89.list, c(train.y89), t, L)

####======================= predict probability ===============================####
## convert matrix into list
# class 0 pair
num <- dim(test.x)[1]
test.x.list <- as.list(1:num)
for (i in 1:num) test.x.list[[i]] <- matrix(test.x[i,], nrow=1,byrow=TRUE)
# class 0 pair
obj201 <- predict.fsvm.prob(obj01, test.x.list)
obj202 <- predict.fsvm.prob(obj02, test.x.list)
obj203 <- predict.fsvm.prob(obj03, test.x.list)
obj204 <- predict.fsvm.prob(obj04, test.x.list)
obj205 <- predict.fsvm.prob(obj05, test.x.list)
obj206 <- predict.fsvm.prob(obj06, test.x.list)
obj207 <- predict.fsvm.prob(obj07, test.x.list)
obj208 <- predict.fsvm.prob(obj08, test.x.list)
obj209 <- predict.fsvm.prob(obj09, test.x.list)
# class 1 pair
obj212 <- predict.fsvm.prob(obj12, test.x.list)
obj213 <- predict.fsvm.prob(obj13, test.x.list)
obj214 <- predict.fsvm.prob(obj14, test.x.list)
obj215 <- predict.fsvm.prob(obj15, test.x.list)
obj216 <- predict.fsvm.prob(obj16, test.x.list)
obj217 <- predict.fsvm.prob(obj17, test.x.list)
obj218 <- predict.fsvm.prob(obj18, test.x.list)
obj219 <- predict.fsvm.prob(obj19, test.x.list)
# class 2 pair
obj223 <- predict.fsvm.prob(obj23, test.x.list)
obj224 <- predict.fsvm.prob(obj24, test.x.list)
obj225 <- predict.fsvm.prob(obj25, test.x.list)
obj226 <- predict.fsvm.prob(obj26, test.x.list)
obj227 <- predict.fsvm.prob(obj27, test.x.list)
obj228 <- predict.fsvm.prob(obj28, test.x.list)
obj229 <- predict.fsvm.prob(obj29, test.x.list)
# class 3 pair
obj234 <- predict.fsvm.prob(obj34, test.x.list)
obj235 <- predict.fsvm.prob(obj35, test.x.list)
obj236 <- predict.fsvm.prob(obj36, test.x.list)
obj237 <- predict.fsvm.prob(obj37, test.x.list)
obj238 <- predict.fsvm.prob(obj38, test.x.list)
obj239 <- predict.fsvm.prob(obj39, test.x.list)
# class 4 pair
obj245 <- predict.fsvm.prob(obj45, test.x.list)
obj246 <- predict.fsvm.prob(obj46, test.x.list)
obj247 <- predict.fsvm.prob(obj47, test.x.list)
obj248 <- predict.fsvm.prob(obj48, test.x.list)
obj249 <- predict.fsvm.prob(obj49, test.x.list)
# class 5 pair
obj256 <- predict.fsvm.prob(obj56, test.x.list)
obj257 <- predict.fsvm.prob(obj57, test.x.list)
obj258 <- predict.fsvm.prob(obj58, test.x.list)
obj259 <- predict.fsvm.prob(obj59, test.x.list)
# class 6 pair
obj267 <- predict.fsvm.prob(obj67, test.x.list)
obj268 <- predict.fsvm.prob(obj68, test.x.list)
obj269 <- predict.fsvm.prob(obj69, test.x.list)
# class 7 pair
obj278 <- predict.fsvm.prob(obj78, test.x.list)
obj279 <- predict.fsvm.prob(obj79, test.x.list)
# class 8 pair
obj289 <- predict.fsvm.prob(obj89, test.x.list)

####============================= Pairwise Calculation =====================================####
###=============================== Algorithm 2. ============================================####
n.test <- length(test.y)
result <- as.list(1:n.test)
K <- 10

# Pairwise Coupling
for(ii in 1:n.test){
  r <- matrix(0, K, K)
  r[lower.tri(r)] <- c(obj201$prob[ii], obj202$prob[ii], obj203$prob[ii], obj204$prob[ii], obj205$prob[ii], obj206$prob[ii], obj207$prob[ii], obj208$prob[ii], obj209$prob[ii], 
      obj212$prob[ii], obj213$prob[ii], obj214$prob[ii], obj215$prob[ii], obj216$prob[ii], obj217$prob[ii], obj218$prob[ii], obj219$prob[ii],
      obj223$prob[ii], obj224$prob[ii], obj225$prob[ii], obj226$prob[ii], obj227$prob[ii], obj228$prob[ii], obj229$prob[ii],
      obj234$prob[ii], obj235$prob[ii], obj236$prob[ii], obj237$prob[ii], obj238$prob[ii], obj239$prob[ii],
      obj245$prob[ii], obj246$prob[ii], obj247$prob[ii], obj248$prob[ii], obj249$prob[ii],
      obj256$prob[ii], obj257$prob[ii], obj258$prob[ii], obj259$prob[ii],
      obj267$prob[ii], obj268$prob[ii], obj269$prob[ii],
      obj278$prob[ii], obj279$prob[ii],
      obj289$prob[ii])
  r[upper.tri(r)] <- t(r)[upper.tri(r)]
  one <- matrix(1,K,K)
  one[lower.tri(one, diag=TRUE)] <- 0
  r <- one - r
  r <- abs(r)
  
  ## Algorithm2.
  ### Create Q matrix
  Q <- matrix(0,K,K)
  for(i in 1:K){
    for(j in 1:K){
      Q[i,j] <- -r[j,i]*r[i,j]
    }
  }
  diag(Q) <- colSums(r^2)
  
  ### (1) Initialize P
  p <- matrix(rep(1/K),K) # 이렇게 초기화해도돼나Q
  p[K] <- 1-sum(p[-K])
  
  ### (2) Repeat (t = 1, 2, 3, ..., K, 1, ...)
  t <- 1
  iter.n <- 1
  
  while(TRUE){
    a <- 1/Q[t,t]
    b <- t(p) %*% Q %*% p
    p[t,] <- a * ( -as.vector(Q[t,-t]) %*% p[-t]  + b)
    
    ## normalize
    p <- p/sum(p)
    
    ## Condition (21) check
    tmp <- Q %*% p
    tmp2 <- matrix(outer(tmp, tmp, "-"), K, K)
    tmp3 <- tmp2[upper.tri(tmp2)]
    # print(c(tmp3,iter.n))
    # print(p)
    idx <- which(max(p) == p)[1]
    c <- c(p[idx] - p[-idx])
    # if (length(unique(sign(tmp))) == 1 && abs(tmp3) < 1e-03)
    #   if(sum(p) == 1 && c > 0.99) break
    if (length(unique(sign(tmp))) == 1 && abs(tmp3) < 1e-05 && sum(p) == 1) break
    
    ## re-indexing t
    if (t == K) {
      t <- 1
    } else {
      t <- (t + 1)
    }
    
    iter.n <- iter.n + 1 # counting the iteration
  }
  
  # Save results
  result[[ii]] <- list(p = p, iter = iter.n, Q = Q)
}

result[[1]]$p
test.y

p.class <- c()
pred.p <- c()

for (ii in 1:n.test){
  p.class[ii] <- which(max(result[[ii]]$p) == result[[ii]]$p)
  pred.p[ii] <- max(result[[ii]]$p)
}

dim(test.y)
length(p.class)
sum(test.y == p.class)
test.y[test.y == p.class]
pred.p





