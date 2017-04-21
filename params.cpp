#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
using namespace Rcpp;
using namespace std;
//[[Rcpp::depends(RcppArmadillo)]]



/*
██   ██ ███████  █████  ██████  ███████ ██████
██   ██ ██      ██   ██ ██   ██ ██      ██   ██
███████ █████   ███████ ██   ██ █████   ██████
██   ██ ██      ██   ██ ██   ██ ██      ██   ██
██   ██ ███████ ██   ██ ██████  ███████ ██   ██
*/

arma::mat log_prob_X(arma::mat X, arma::mat D, arma::mat P );
arma::mat log_prob_D(arma::mat D, arma::mat M, arma::vec phi);


/* ------------------------ basic ----------------------------------*/
inline bool any_cpp(LogicalVector lv)
{
  return is_true(any(lv));
}

arma::mat matrix_sub(arma::mat M, LogicalVector a, int b)
{
  // b=1: select row
  // b=2: select column
  arma::mat out;
  if(b==2){
    arma::colvec z=as<arma::colvec>(a);
    out=M.cols(find(z==1));
  } else if(b==1){
    arma::rowvec z=as<arma::rowvec>(a);
    out=M.rows(find(z==1));
  }

  return out;
}


/*
████████ ██████  ███████ ███████
   ██    ██   ██ ██      ██
   ██    ██████  █████   █████
   ██    ██   ██ ██      ██
   ██    ██   ██ ███████ ███████
*/


//----------------------------------------------------------------------------------------
/*  -------------------------------- tree related  -------------------------------*/
//----------------------------------------------------------------------------------------

//[[Rcpp::export]]
LogicalVector find_desc_cpp(int k, NumericVector tree)
{
  // returning 1,0 vector indicating whether each subclone is descendant
  LogicalVector parent_set(tree.size());

  parent_set[k-1]=1;
  if (k<tree.size()){
    for (int i=k+1;i<=tree.size();i++){
      parent_set[i-1] = (parent_set[tree[i-1]-1]) ? 1:0;
    }
  }
  return parent_set;
}



//[[Rcpp::export]]
LogicalVector find_child_cpp(double k, NumericVector tree)
{
  return tree==k;
}


//[[Rcpp::export]]
NumericMatrix Lo_to_L_cpp(NumericMatrix B, NumericVector tree)
{
  int J=B.nrow();// number of muts
  int K=tree.length(); // number of subclones

  NumericMatrix out(J,K);
  NumericVector a(2);
  NumericVector b(K);

  LogicalVector temp(K);
  for (int i=0;i<J;i++){
    a=B.row(i);
    b=wrap(2*arma::ones(K));

    if(a[0]!=0){
      temp=find_desc_cpp(a[0], tree);
      b[temp] =  a[1]+2 ;
    }

    out.row(i)=b;
  }
  return out;
}

//[[Rcpp::export]]
NumericMatrix L_to_Lo_cpp(NumericMatrix B)
{
  int J=B.nrow();// number of muts
  int K=B.ncol(); // number of subclones

  NumericMatrix out(J,2);

  NumericVector temp(2);
  for (int j=0;j<J;j++){
    temp[0]=temp[1]=0;
      for(int k=0;k<K;k++){
        if(B(j,k) != 2){temp[1]= B(j,k)-2;temp[0]=k+1;break;}
      }
    out.row(j)=temp;
  }
  return out;
}


//[[Rcpp::export]]
NumericMatrix Zo_to_Z_cpp(NumericMatrix B, NumericVector tree)
{
  int J=B.nrow();// number of muts
  int K=tree.length(); // number of subclones

  NumericMatrix out(J,K);
  NumericVector a(2);
  NumericVector b(K);

  LogicalVector temp(K);
  for (int i=0;i<J;i++){
    a=B.row(i);
    b=wrap(arma::zeros(K));
    temp=find_desc_cpp(a[0], tree);
    b[temp] =  a[1] ;
    out.row(i)=b;
  }
  return out;
}

//[[Rcpp::export]]
NumericMatrix Z_to_Zo_cpp(NumericMatrix B)
{
  int J=B.nrow();// number of muts
  int K=B.ncol(); // number of subclones

  NumericMatrix out(J,2);

  NumericVector temp(2);
  for (int j=0;j<J;j++){
    temp[0]=temp[1]=0;
    for(int k=0;k<K;k++){
      if(B(j,k) != 0){temp[1]= B(j,k);temp[0]=k+1;break;}
    }
    out.row(j)=temp;
  }
  return out;
}


/*
███████  █████  ███    ███ ██████  ██      ███████     ███████
██      ██   ██ ████  ████ ██   ██ ██      ██          ██
███████ ███████ ██ ████ ██ ██████  ██      █████       █████
     ██ ██   ██ ██  ██  ██ ██      ██      ██          ██
███████ ██   ██ ██      ██ ██      ███████ ███████     ██
*/



//------------------------------------------------------------------------
// ------------------------- sampling of theta --------------------------
//------------------------------------------------------------------------

//[[Rcpp::export]]
double prior_theta(arma::vec theta_t, List Params)
{
  double r=Params["r"];
  double out= (r-1)*sum(log(theta_t)) - sum(theta_t);
  return out;
}

//[[Rcpp::export]]
double log_p_theta_cpp(arma::vec theta_t,double phi_t, List samp, arma::vec xt,
                       arma::vec dt, List Params, double temper)
{
  double out;

  arma::mat Z=as<arma::mat>(samp["Z"]), L=as<arma::mat>(samp["L"]);
  // double overdispersion = samp["over_omega"];
  double Gt = sum(theta_t);
  arma::mat m = L*(theta_t/Gt);
  arma::mat p = Z*(theta_t/Gt)/m;

  arma::vec phi(1);
  phi[0] = phi_t;
  arma::mat px = log_prob_X(xt,dt,p);
  arma::mat pd = log_prob_D(dt,m,phi);

  out =  accu(px) + accu(pd) + prior_theta(theta_t, Params);

  return out/temper;
}




NumericVector samp_theta_cpp(arma::vec theta_t, int i, double phi_t, List samp,
                    arma::vec xt, arma::vec dt, List Params,double p0, double tune_par,
                    double temper, double err)
{ 
  NumericVector out(2);

  double theta0= theta_t[i];

  //propose
  double theta1= R::rgamma(tune_par*theta0, 1/tune_par);
  if(theta1==0){theta1 += err;}



  double p_01= R::dgamma(theta1,tune_par*theta0,1/tune_par,0);
  double p_10= R::dgamma(theta0,tune_par*theta1,1/tune_par,0);


  //posterior of new data

  if(p0==0){p0=log_p_theta_cpp(theta_t,phi_t,samp,xt,dt,Params,temper);}
  arma::vec new_theta_t=theta_t;
  new_theta_t[i]=theta1;

  double p1=log_p_theta_cpp(new_theta_t,phi_t,samp,xt,dt,Params,temper);
  // acceptance probability

  double temp_prob=exp(p1-p0)*(p_10/p_01);
  double acc_prob= (temp_prob >= 1)? 1 : temp_prob;


  double u=R::runif(0,1);
  if(u<=acc_prob){
    out[0]=theta1;
    out[1]=p1;
  } else{
    out[0]=theta0;
    out[1]=p0;
  }

  return out;
}


//[[Rcpp::export]]
arma::mat samp_theta_all (List samp, List Params, arma::mat X, arma::mat D,
                             double theta_tune,double temper)
{  
  arma::mat out=as<arma::mat>(samp["theta"]);
  int K=out.n_rows, T0=out.n_cols;
  arma::vec phi= as<arma::vec>(Params["phi"]);

  double p0=0;
  double phi_t;


  for (int j =0; j<T0; j++){
    phi_t=phi.at(j);
    for (int i=0; i<K;i++){
      NumericVector update= samp_theta_cpp(out.col(j),i,phi_t,samp,X.col(j),D.col(j),
                                Params,p0,theta_tune,temper,0.01);
      out(i,j)=update[0];
      p0=update[1];
      // update vector
      if(i==(K-1)) {p0=0;}
    }

  }

  return out;
}

/*
███████  █████  ███    ███ ██████  ██      ███████     ██
██      ██   ██ ████  ████ ██   ██ ██      ██          ██
███████ ███████ ██ ████ ██ ██████  ██      █████       ██
     ██ ██   ██ ██  ██  ██ ██      ██      ██          ██
███████ ██   ██ ██      ██ ██      ███████ ███████     ███████
*/



//----------------------------------------------------------------
// ------------------------- sampling of L--------------------------
//--------------------------------------------------------------------

double prior_L_cpp(NumericVector Lj, List samp, List Params)
{
  double out;

  double K=Params["K"];
  double pi=samp["pi"];
  double max_CN=Params["max_CN"];

  LogicalVector temp= (Lj != 2);
  if(any_cpp(temp)){
    out=(1-pi)/((K-1)*max_CN);
  } else {
    out=pi;
  }

  return log(out);
}


//[[Rcpp::export]]
double log_p_L_cpp(arma::vec Lj, int j, arma::mat D, arma::mat X, List samp, List Params, double temper )
{
  double out=0;

  arma::mat Frac = as<arma::mat>(samp["Frac"]);
  arma::mat Z =  as<arma::mat>(samp["Z"]);
  arma::vec phi = as<arma::vec>(Params["phi"]);
  
  arma::mat ave_CN = (Lj.t() * Frac);  // row vector of average CN
  arma::mat ave_mut = (Z.row(j-1) * Frac) ; // row vector of average mut
  arma::mat p = ave_mut/ave_CN;

  arma::mat px = log_prob_X(X.row(j-1),D.row(j-1),p);
  arma::mat pd = log_prob_D(D.row(j-1),ave_CN,phi);

  out = accu(pd) + accu(px);
  out /= temper;
  return out;
}

//[[Rcpp::export]]
NumericVector propose_L_cpp(NumericVector tree, arma::mat Z_seg, List Params)
{
  int K=Params["K"];
  int max_CN=Params["max_CN"];

  int origin=R::runif(2,K+1); // randomly select an origin
  LogicalVector desc=find_desc_cpp(origin, tree);

  arma::mat temp=matrix_sub(Z_seg,desc,2);
  int min_cnv=temp.max(); // make sure L>Z
  int new_cn= R::runif(min_cnv,max_CN+1); // new copy number

  NumericVector out(K,2.0);
  out[desc]=new_cn;

  return out;
}


//[[Rcpp::export]]
NumericVector Gibbs_L_cpp(arma::uvec indVec, List samp, arma::mat X, arma::mat D,
                               List Params, double temper)
{
  arma::mat L=samp["L"], Z=samp["Z"];
  int K= Params["K"];
  int maxCN=Params["max_CN"];
  NumericVector tree=samp["Ttree"];
  arma::mat Z_seg=Z.rows(indVec-1);


  arma::vec p_vec=arma::zeros<arma::vec>((K-1)*maxCN+1);
  double p0;

  // k=0
  int min_cnv=Z_seg.max(); // make sure L>Z
  arma::vec normal_l(K);
  normal_l.fill(2);

  if(min_cnv<=2){
    p0= prior_L_cpp(wrap(normal_l),samp,Params);
    for(arma::uvec::iterator it_j=indVec.begin(); it_j != indVec.end(); it_j++){
      p0 += log_p_L_cpp(normal_l,*it_j,D,X,samp,Params,temper);
    }
    p_vec.at(0)=p0;
  }

  
  int ind=1;
  arma::vec new_l(K);
  for (int k=2; k<=K; k++){
    for(int c=-2; c<=maxCN-2 ; c++){

      if(c==0) continue;
      LogicalVector desc=find_desc_cpp(k, tree);
      min_cnv=matrix_sub(Z_seg,desc,2).max();
       if(c+2>=min_cnv){
         // create l vector
         new_l.fill(2);
         new_l.elem(find(as<arma::uvec>(desc))) = new_l.elem(find(as<arma::uvec>(desc))) +c;
         // calculate posterior
         p0=prior_L_cpp(wrap(new_l),samp,Params);
         for(arma::uvec::iterator it_j=indVec.begin(); it_j != indVec.end(); it_j++){
           p0 += log_p_L_cpp(new_l,*it_j,D,X,samp,Params,temper);
         }
         p_vec.at(ind)=p0;
       }
       ind++;
    }
  }

  double increase= p_vec.elem(find(p_vec)).max();
  p_vec.elem(find(p_vec))= exp( p_vec.elem(find(p_vec)) - increase );

  IntegerVector frame=Range(0,(K-1)*maxCN);
  int ind_samp= as<int>(Rcpp::RcppArmadillo::sample(frame,1,1,p_vec));
  if (ind_samp==0) {return wrap(normal_l);}


  int k=2+(ind_samp-1)/maxCN;
  int c=(ind_samp-1)-(k-2)*maxCN-2;
  c= (c>=0) ? c+1 : c;

  arma::vec new_lo(2);
  new_lo.at(0)= k;
  new_lo.at(1)= c;

  NumericVector out = Lo_to_L_cpp(wrap(new_lo.t()),tree);
  return out;
}


//[[Rcpp::export]]
NumericMatrix samp_L_all(List samp, arma::mat X, arma::mat D, List Params, double temper)
{
  int J=Params["J"]; // number locus
  int K=Params["K"]; // number of subclones
  arma::mat segs=Params["segments"];
  int nsegs= segs.n_rows;

  NumericMatrix out(J,K);
  NumericVector temp(K);


  for (int i=0; i<nsegs; i++){
    // generate sequence for this segment
    int start=segs(i,0), end=segs(i,1);
    arma::uvec one_seg(end-start+1);
    for(arma::uvec::iterator uvec_it=one_seg.begin(); uvec_it != one_seg.end(); uvec_it++){
      *uvec_it=start++;
    }

    temp=Gibbs_L_cpp(one_seg,samp,X,D,Params,temper);

    // assign output to all rows in this segment
    for(arma::uvec::iterator uvec_it=one_seg.begin(); uvec_it != one_seg.end(); uvec_it++){
    out.row(*uvec_it-1)=temp;
    }
  }

  return out;
}


/*
███████  █████  ███    ███ ██████  ██      ███████     ███████
██      ██   ██ ████  ████ ██   ██ ██      ██             ███
███████ ███████ ██ ████ ██ ██████  ██      █████         ███
     ██ ██   ██ ██  ██  ██ ██      ██      ██           ███
███████ ██   ██ ██      ██ ██      ███████ ███████     ███████
*/



//---------------------------------------------------------------------------
// ------------------------- sampling of Z  --------------------------
//------------------------------------------------------------------------

double prior_Z_cpp(NumericVector Zj, List Params)
{
  int a=max(Zj);
  int max_mut=Params["max_mut"];
  double zeta=Params["zeta"];

  double p_a=pow(zeta,a);
  double p_all=(1-pow(zeta,max_mut+1))/(1-zeta);
  p_all -= 1;

  double out=p_a/p_all;
  return log(out);
}


//[[Rcpp::export]]
double log_p_Z_cpp(arma::vec Zj,int j,  arma::mat D, arma::mat X, List samp, List Params,
                   double temper)
{
  arma::mat M=samp["M"];
  arma::mat Frac=samp["Frac"];
  arma::mat ave_mut= Zj.t() * Frac; // row vector for average mutation
  arma::mat p = ave_mut / M.row(j-1);

  arma::mat px = log_prob_X(X.row(j-1),D.row(j-1),p);

  double out = accu(px) + prior_Z_cpp(wrap(Zj),Params);

  out = out/temper;
  return out;
}


//[[Rcpp::export]]
NumericVector propose_Z_cpp(NumericVector tree, NumericVector Lj, List Params)
{
  int K=Params["K"];
  // sample an origin with L > 0
  LogicalVector loc=(Lj > 0);
  loc[0]=0; // change value for normal clone

  int s=sum(loc);
  int u=R::runif(1,s+1);
  int cumsum=0, index=1;
  for (LogicalVector::iterator it=loc.begin(); it !=loc.end(); it++){
    cumsum += loc[index-1];
    if(cumsum==u)break;
    index++;
  }


  LogicalVector desc=find_desc_cpp(index,tree); // subclones to be changed

  //sample number of mutant copy
  int max_mut=Params["max_mut"];
  NumericVector Lj_desc=Lj[desc];
  int min_Lj_desc=min(Lj_desc);

  int max_mut2 = (max_mut > min_Lj_desc) ? min_Lj_desc : max_mut; // maximum mutant copy allowed

  u=R::runif(1,max_mut2+1);

  NumericVector out(K,0.0);
  out[desc]=u;
  return out;
}



//[[Rcpp::export]]
NumericVector Gibbs_Z_cpp(int j, arma::vec Lj, List samp, arma::mat X, arma::mat D,
                          List Params, double temper)
{
  NumericVector tree=samp["Ttree"];

  int K=Params["K"];
  int maxMut=Params["max_mut"];

  arma::vec p_vec=arma::zeros<arma::vec>((K-1)*maxMut);
  int ind=0;
  for(int k=2; k<=K; k++){
    for(int c=1; c<=maxMut; c++){
      //prob for (k,c)
      arma::uvec q=find(as<arma::uvec>(find_desc_cpp(k,tree)));
      int temp_max=Lj.elem(q).min();
      if(temp_max<c) {ind++;continue;}

      arma::vec cur_z0(2);
      cur_z0.at(0)=k;
      cur_z0.at(1)=c;

      NumericMatrix temp_Z=Zo_to_Z_cpp(wrap(cur_z0.t()),tree);
      p_vec.at(ind)= log_p_Z_cpp(as<arma::vec>(temp_Z),j,D,X,samp,Params,temper);
      ind++;
    }
  }

  double increase= p_vec.elem(find(p_vec)).max();
  p_vec.elem(find(p_vec))= exp( p_vec.elem(find(p_vec)) - increase );


  IntegerVector frame=Range(1,(K-1)*maxMut);
  int ind_samp= as<int>(Rcpp::RcppArmadillo::sample(frame,1,1,p_vec));
  arma::vec new_zo(2);
  int k=2+(ind_samp-1)/maxMut;
  int c=ind_samp-(k-2)*maxMut;
  new_zo.at(0)= k;
  new_zo.at(1)= c;

  NumericVector new_z = Zo_to_Z_cpp(wrap(new_zo.t()),tree);
  return new_z;
}



//[[Rcpp::export]]
NumericMatrix samp_Z_all(List samp, arma::mat X, arma::mat D, List Params, double temper)
{
  int J=Params["J"], K=Params["K"];
  arma::mat Z=samp["Z"], L=samp["L"];
  NumericMatrix out(J,K);
  NumericVector temp(K);

  arma::vec Zj,Lj;
  for (int j=1; j<=J; j++){
    Zj=Z.row(j-1).t();
    Lj=L.row(j-1).t();
    temp=Gibbs_Z_cpp(j,Lj,samp,X,D,Params,temper);
    out.row(j-1)=temp;
  }
  return out;
}


/*
 ██████  ████████ ██   ██ ███████ ██████  ███████
██    ██    ██    ██   ██ ██      ██   ██ ██
██    ██    ██    ███████ █████   ██████  ███████
██    ██    ██    ██   ██ ██      ██   ██      ██
 ██████     ██    ██   ██ ███████ ██   ██ ███████
*/


// ------------------------------------------------------------
/*  --------------------- others ---------------------- */
// ------------------------------------------------------------


//[[Rcpp::export]]
arma::mat log_prob_X(arma::mat X, arma::mat D, arma::mat P )
{
  int J = X.n_rows;
  int t = X.n_cols;
  arma::mat out(J,t );
  for (int i=0; i<J; i++){
    for(int j=0; j<t; j++){
      out.at(i,j) = R::dbinom(X.at(i,j),D.at(i,j),P.at(i,j),1);
    }
  }
  return out;
}

//[[Rcpp:export]]
arma::mat log_prob_D(arma::mat D, arma::mat M, arma::vec phi)
{
  int J = D.n_rows;
  int t = D.n_cols;
  arma::mat out(J,t);
  for (int i=0; i<J; i++){
    for(int j=0; j<t; j++){
      double lambda=phi[j]*M.at(i,j)/2;
      out.at(i,j)= R::dpois(D.at(i,j),lambda,1);
    }
  }
  return out;
}


// calculate prior for all parameters
//[[Rcpp::export]]
double log_prior_all(List samp, List Params, double temper)
{
  int T0=Params["T0"];
  int J=Params["J"];
  arma::mat L=samp["L"], Z=samp["Z"], theta=samp["theta"];
  // double omega = samp["over_omega"];

  int a=Params["a_pi"];
  int b=Params["b_pi"];
  double pi=samp["pi"];

  double out=0;
// pi
  double temp=R::dbeta(pi,a,b,1);
  out += temp;

//theta prior
  for (int t=0; t<T0; t++){
    out += prior_theta( theta.col(t) ,Params);
  }

//L prior
  for(int j=0; j< L.n_rows; j++){
    out += prior_L_cpp( wrap(L.row(j).t()) ,samp,Params);
  }

//Z prior
  for(int j=0; j<J; j++){
    out += prior_Z_cpp(wrap(Z.row(j).t()),Params);
  }

// omega prior
  // out +=  log_prior_omega(omega, Params);
  
  return out/temper;
}


// calculate likilihood matrix l_jt=p(x_jt | ...) p(d_jt | ...)
//[[Rcpp::export]]
arma::mat likelihood_mat(arma::mat X, arma::mat D, arma::mat Frac, arma::mat L, arma::mat Z, double over_dispersion, List Params)
{
  arma::vec phi=Params["phi"];
  arma::mat M = L * Frac;

  arma::mat P= Z * Frac / M;
  arma::mat out = log_prob_X(X,D,P) + log_prob_D(D,M,phi);

  return out;
}
