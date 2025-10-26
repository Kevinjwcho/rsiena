# install.packages("devtools")
library(devtools)

# 패키지 로드
load_all()   # 현재 디렉토리의 R/ 아래 파일들 불러오기

set.seed(7)
X <- array(sample(c(0,1,10,11,NA), 8*8*8*3, replace=TRUE, prob=c(.86,.09,.02,.02,.01)),
           dim=c(8,8,8,3))
Y <- sienaDependent(X, type = "threeway")
attr(Y, "type")      # "threeway"
attr(Y, "netdims")   # 8 8 8 3
attr(Y, "nodeSet")   # "Actors"
