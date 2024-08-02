test_that("iccmulti() estimates binary ICC", {
  set.seed(1234)
  dat2 = rccat(rho=0.2,prop=0.3,noc=20,csize=20)
  expect_equal(round(as.numeric(iccmulti(cid,y,dat2,binmethod="rm",ci.type="rm")$estimates[2]),3),
               round(as.numeric(ICCbin::iccbin(cid,y,dat2,method="rm",ci.type="rm")$estimates[2]),3)
               )
})
