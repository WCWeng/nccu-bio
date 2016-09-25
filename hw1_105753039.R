# read PAM1 from data
temp <- read.table("pam1.txt")
pam1 <- as.matrix(temp)

# check PAM1 data
dim(pam1)
str(pam1)

# construct PAM250 from PAM1
pam1 <- pam1 / 10000.0 #轉為機率
pam250 <- pam1
for (x in c(1:249)) pam250 <- pam250 %*% pam1
pam250 <- pam250 * 100.0 #轉為百分比機率

# output PAM250 as a file
write.table(pam250, "pam250.txt")
