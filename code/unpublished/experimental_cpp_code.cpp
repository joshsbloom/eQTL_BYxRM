calc.mr2.code.noInv <- '   arma::mat XYr = Rcpp::as<arma::mat>(XYr_);
   int n = XYr.n_rows; int p = XYr.n_cols;
   arma::mat XXr = Rcpp::as<arma::mat>(XXr_);
   arma::mat r2a(p,3, arma::fill::zeros); 
   int mspace = as<int>(mspace_);
   int ii;
   int jj;
   int kk;
   double r2;
   double det2;
   
   arma::uvec i2 = arma::zeros<arma::uvec>(2);
   arma::mat RRinvA(2,2, arma::fill::zeros);
   arma::mat RRinvAtmp(2,2, arma::fill::zeros);

   for(ii=0; ii<(n-mspace); ii++) {
      //Rcpp::Rcout << ii <<  std::endl;
      for(jj=(ii+mspace); jj<n; jj++){
          i2(0)=ii; i2(1)=jj;
          arma::mat X = XYr.rows(i2);
          RRinvAtmp =XXr.submat(i2,i2);
          det2 = 1/((RRinvAtmp(0,0)*RRinvAtmp(1,1)) - (RRinvAtmp(0,1)*RRinvAtmp(1,0)));
          RRinvA(0,0)=RRinvAtmp(1,1)*det2;
          RRinvA(0,1)=-RRinvAtmp(0,1)*det2;
          RRinvA(1,0)=-RRinvAtmp(1,0)*det2;
          RRinvA(1,1)=RRinvAtmp(0,0)*det2;

          //arma::mat RRinv = inv_sympd(XXr.submat(i2,i2));
          
          for(kk=0; kk<p; kk++) {
              r2 = arma::as_scalar(trans(X.col(kk)) * RRinvA * X.col(kk));
              if(r2>r2a(kk,2)) {
                r2a(kk,0)=ii+1;
                r2a(kk,1)=jj+1;
                r2a(kk,2)=r2;
              }
           }
    }
    }    
   return Rcpp::wrap(r2a);
'
scantwo.additive.noInv <- cxxfunction(signature(XYr_="matrix", XXr_="matrix", mspace_="numeric"),
                                    calc.mr2.code.noInv, plugin="RcppArmadillo")


calc.mr3.code.noInv=
'
arma::mat XYr = Rcpp::as<arma::mat>(XYr_);
   int n = XYr.n_rows; int p = XYr.n_cols;
   arma::mat XXr = Rcpp::as<arma::mat>(XXr_);
   arma::mat r2a(p,4, arma::fill::zeros); 
   int mspace = as<int>(mspace_);
   int ii;
   int jj;
   int ll;
   int kk;
   double r2;
   double u;
   double det2;
   arma::uvec i2 = arma::zeros<arma::uvec>(2);
   arma::uvec i3 = arma::zeros<arma::uvec>(3);
   arma::uvec i4 = arma::zeros<arma::uvec>(1);
   arma::mat RRinv(3,3, arma::fill::zeros);
   arma::mat RRinvA(2,2, arma::fill::zeros);
   arma::mat RRinvAtmp(2,2, arma::fill::zeros);
   for(ii=0; ii<(n-mspace-mspace); ii++) {
      //Rcpp::Rcout << ii <<  std::endl;
      for(jj=(ii+mspace); jj<(n-mspace); jj++){
         i2(0)=ii; i2(1)=jj;
         i3(0)=ii; i3(1)=jj; 
         
         RRinvAtmp =XXr.submat(i2,i2);
         det2 = 1/((RRinvAtmp(0,0)*RRinvAtmp(1,1)) - (RRinvAtmp(0,1)*RRinvAtmp(1,0)));
         RRinvA(0,0)=RRinvAtmp(1,1)*det2;
         RRinvA(0,1)=-RRinvAtmp(0,1)*det2;
         RRinvA(1,0)=-RRinvAtmp(1,0)*det2;
         RRinvA(1,1)=RRinvAtmp(0,0)*det2;

         //RRinvA = inv_sympd(XXr.submat(i2,i2));
          for(ll=(jj+mspace); ll<n; ll++){
              //Rcpp::Rcout << ll <<  std::endl;
              i3(2)=ll;
              i4(0)=ll;
              arma::vec delta = trans(XXr.submat(i4,i2));
              arma::mat idm  = RRinvA * delta;
              u = arma::as_scalar(1 - trans(delta) * idm);
              RRinv(arma::span(0,1), arma::span(0,1)) = RRinvA + (idm * trans(delta) * RRinvA)/u;
              RRinv(2, arma::span(0,1)) = trans(-idm/u) ;
              RRinv(arma::span(0,1),2) = -idm/u ;
              RRinv(2,2)=1/u;
              
              arma::mat X = XYr.rows(i3);
              //Rcpp::Rcout << RRinv <<  std::endl;
              for(kk=0; kk<p; kk++) {
                  r2 = arma::as_scalar(trans(X.col(kk)) * RRinv * X.col(kk));
                  //Rcpp::Rcout << r2 <<  std::endl;

                  if(r2>r2a(kk,3)) {
                    r2a(kk,0)=ii+1;
                    r2a(kk,1)=jj+1;
                    r2a(kk,2)=ll+1;
                    r2a(kk,3)=r2;
                  }
               }
               
        }
    }
    }    
   return Rcpp::wrap(r2a);
'
scanthree.additive.noInv <- cxxfunction(signature(XYr_="matrix", XXr_="matrix", mspace_="numeric"),
                                                  calc.mr3.code.noInv, plugin="RcppArmadillo")



colMaxIndRcpp <- cxxfunction(signature(X_="numeric"), plugin="Rcpp", body='
   Rcpp::NumericMatrix X(X_);
   int n = X.ncol();
   Rcpp::NumericVector V(n);
   for (int i=0; i<n; i++) {
      Rcpp::NumericVector W = X.column(i);
      Rcpp::NumericVector::iterator it = std::max_element(W.begin(), W.end());
      V[i] = it - W.begin() +1;  
   }
   return(V);
 ')

colMaxRcpp <- cxxfunction(signature(X_="numeric"), plugin="Rcpp", body='
   Rcpp::NumericMatrix X(X_);
   int n = X.ncol();
   Rcpp::NumericVector V(n);
   for (int i=0; i<n; i++) {
      Rcpp::NumericVector W = X.column(i);
      V[i] = *std::max_element(W.begin(), W.end());  // from the STL
   }
   return(V);
 ')



#svXXr=svd(XXr)
#svtXXr=svd(stpreds%*%t(stpreds))
#lmpca=list()
#for(i in 1:20) {
#    print(i)
#    print(summary(lm(sy[,1]~svtXXr$u[,i]-1)))
#}
#testpcr=pcr(sy[,1]~stpreds, 1)
calc.mr2.code.r2only <- '
   arma::mat X = Rcpp::as<arma::mat>(X_);
   int n = X.n_rows; int p = X.n_cols;
   arma::mat rinv = Rcpp::as<arma::mat>(rinv_);
   arma::vec r2 = Rcpp::as<arma::vec>(r2_);
   int ii;
   for(ii=0; ii<p; ii++) {
       arma::mat cp1 = rinv * X.col(ii);
       arma::mat cp2 = trans(X.col(ii)) * cp1;
      r2[ii] = cp2[0,0];
   }
   return Rcpp::wrap(r2);
  '
mr2.r2 <- cxxfunction(signature(X_ ="matrix", rinv_="matrix", r2_="numeric"),
                                                  calc.mr2.code.r2only, plugin="RcppArmadillo")

# n is markers
# p is phenotypes
#//int nmar = as<int>(nmar_);
#//          arma::mat rinv = inv_sympd(RR);
#//Rcpp::Rcout << ii <<  std::endl;
#//Rcpp::Rcout << ii <<  std::endl;
#//arma::uvec isgreater= arma::find(r2>r2a.col(2));
#          //r2a.rows(isgreater).col(0) = r2.elem(isgreater);
calc.mr2.code <- '
   arma::mat XYr = Rcpp::as<arma::mat>(XYr_);
   int n = XYr.n_rows; int p = XYr.n_cols;
   arma::mat XXr = Rcpp::as<arma::mat>(XXr_);
   arma::mat r2a(p,3, arma::fill::zeros); 
   int mspace = as<int>(mspace_);
   int ii;
   int jj;
   int kk;
   double r2;
   arma::uvec i2 = arma::zeros<arma::uvec>(2);
   for(ii=0; ii<(n-mspace); ii++) {
      //Rcpp::Rcout << ii <<  std::endl;
      for(jj=(ii+mspace); jj<n; jj++){
          i2(0)=ii; i2(1)=jj;
          arma::mat X = XYr.rows(i2);
          arma::mat RRinv = inv_sympd(XXr.submat(i2,i2));
          
          for(kk=0; kk<p; kk++) {
              r2 = arma::as_scalar(trans(X.col(kk)) * RRinv * X.col(kk));
              if(r2>r2a(kk,2)) {
                r2a(kk,0)=ii+1;
                r2a(kk,1)=jj+1;
                r2a(kk,2)=r2;
              }
           }
    }
    }    
   return Rcpp::wrap(r2a);
'
scantwo.additive <- cxxfunction(signature(XYr_="matrix", XXr_="matrix", mspace_="numeric"),
                                    calc.mr2.code, plugin="RcppArmadillo")



calc.mr3.code <- ' arma::mat XYr = Rcpp::as<arma::mat>(XYr_);
   int n = XYr.n_rows; int p = XYr.n_cols;
   arma::mat XXr = Rcpp::as<arma::mat>(XXr_);
   arma::mat r2a(p,4, arma::fill::zeros); 
   int mspace = as<int>(mspace_);
   int ii;
   int jj;
   int ll;
   int kk;
   double r2;
   arma::uvec i2 = arma::zeros<arma::uvec>(3);
   for(ii=0; ii<(n-mspace-mspace); ii++) {
      for(jj=(ii+mspace); jj<(n-mspace); jj++){
        for(ll=(jj+mspace); ll<n; ll++){
          i2(0)=ii; i2(1)=jj;i2(2)=ll;
          arma::mat X = XYr.rows(i2);
          arma::mat RRinv = inv_sympd(XXr.submat(i2,i2));
          for(kk=0; kk<p; kk++) {
              r2 = arma::as_scalar(trans(X.col(kk)) * RRinv * X.col(kk));
              if(r2>r2a(kk,3)) {
                r2a(kk,0)=ii+1;
                r2a(kk,1)=jj+1;
                r2a(kk,2)=ll+1;
                r2a(kk,3)=r2;
              }
           }
        }
    }
    }    
   return Rcpp::wrap(r2a); '
scanthree.additive <- cxxfunction(signature(XYr_="matrix", XXr_="matrix", mspace_="numeric"),
                                                  calc.mr3.code, plugin="RcppArmadillo")

 
#mr3.scantwo.BI.additive <- cxxfunction(signature(XYr_="matrix", XXr_="matrix", mspace_="numeric"),
#                                                  calc.mr3.BI.code, plugin="RcppArmadillo")


