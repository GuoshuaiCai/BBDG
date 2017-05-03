dyn.load("BBDG_MLE.so")

bb.mle <- function(x,which=NULL,method="CG",coe=NULL,group=NULL,pi.m=NULL,ref=NULL,len.r=NULL,a0=a0,afix=afix,d.index=NULL) {

	if(which=="mu"){
        
		out <- optimize(f = .opl.u.i.c,
			lower=0,
			upper =1,
			coe=coe,
			P = a0,
			x = x
			)

		return(list(mu=out$minimum, fval=-out$objective))
	}

	if(which=="d"){
        
		out <- optimize(f = .opl.u.i.c.d,
			lower=0,
			upper =10,
			P = a0,
			x = x
			)
		return(list(d=out$minimum, fval=-out$objective))
	}
	
	if(which=="D_SAS"){
		out <- optim(par = a0,
                     	fn = .opl.u.c.c,
			gr = .opl.g.c.c,
                     	method = "L-BFGS-B",
			control = list(trace = 0,REPORT = 1,maxit = 2000,factr = 1e-22,pgtol = 1e-22),
			lower = c(1e-5,1e-5),
                     	upper = c(100,100),
			index=d.index,
			P=afix,
                     	x = x)
		return(list(d=out$par, fval=-out$value, code=out$message))
	}

	if(which=="theta"){
		if(method=="CG"){
			out <- optim(par = a0,
				     fn = .opl.u.c.D,
				     gr = .opl.g.c.D,
				     method = "CG",
				     control = list(trace = 0,
						    REPORT = 1,
						    maxit = 2000,
						    factr = 1e-22,
						    pgtol = 1e-22),
				     x = x)
		}
		if(method=="L-BFGS-B"){
			out <- optim(par = a0,
				     fn = .opl.u.c.D,
				     gr = .opl.g.c.D,
				     method = "L-BFGS-B",
				     control = list(trace = 0,REPORT = 1,maxit = 2000,factr = 1e-22,pgtol = 1e-22),
				     lower = c(rep(-2, length(a0)-2),0,0),
				     upper = c(rep(2, length(a0)-2),1,20),
				     x = x) 
		}
		return(list(mt=out$par, fval=-out$value, code=out$message))
	}

	if(which=="onepar"){
			out <- optimize(f = .opl.u.c.One,
					lower=0,
					upper =20,
					P = a0,					
					x = x
					)
			return(list(theta=out$minimum, fval=-out$objective))
	}

	if(which=="theta.coe.fix"){
		if(method=="CG"){
			out <- optim(par = a0,
				     fn = .opl.u.coe.fix.c,
				     gr = .opl.g.coe.fix.c,
				     method = "CG",
				     control = list(trace = 0,
						    REPORT = 1,
						    maxit = 2000,
						    factr = 1e-22,
						    pgtol = 1e-22),
				     coe=coe,
				     group=group,
				     pi.m=pi.m,
				     ref=ref,
				     len.r=len.r,
				     x = x)
		}
		if(method=="L-BFGS-B"){
			out <- optim(par = a0,
				     fn = .opl.u.coe.fix.c,
				     gr = .opl.g.coe.fix.c,
				     method = "L-BFGS-B",
				     control = list(trace = 0,REPORT = 1,maxit = 2000,factr = 1e-22,pgtol = 1e-22),
				     lower = c(1e-5,0,0),
				     upper = c(10000,1,100),
				     coe=coe,
				     group=group,
				     pi.m=pi.m,
				     ref=ref,
				     len.r=len.r,
				     x = x)
		}
		return(list(mt=out$par, fval=-out$value, code=out$message))
	}
}


funcv.c<-function(aspect, s, mu, xi, ni, c,d,p,o,index=NULL){

	n<-length(s)
	if(aspect=="g"){
		cval<-.C("bbmle_c_g",as.double(s),as.integer(n),as.double(mu),as.integer(xi),as.integer(ni),as.double(c),as.double(d),as.double(p),as.double(o),result = double(3))
	}
	if(aspect=="u"){
		cval<-.C("bbmle_c_u",as.double(s),as.integer(n),as.double(mu),as.integer(xi),as.integer(ni),as.double(c),as.double(d),as.double(p),as.double(o),result = double(1))
	}
	if(aspect=="i"){
		cval<-.C("bbmle_c_i",as.double(s),as.integer(n),as.double(mu),as.integer(xi),as.integer(ni),as.double(c),as.double(d),as.double(p),as.double(o),result = double(1))
	}
	if(aspect=="du"){
		cval<-.C("bbmle_c_du",as.double(s),as.integer(n),as.double(mu),as.integer(xi),as.integer(ni),as.double(c),as.double(d),as.double(p),as.double(o),result = double(1))
	}
	if(aspect=="dg"){
		cval<-.C("bbmle_c_dg",as.double(s),as.integer(n),as.double(mu),as.integer(xi),as.integer(ni),as.double(c),as.double(d),as.double(p),as.double(o),as.integer(index),result = double(2))
	}

	cval[["result"]]

}

.opl.u.i.c <- function(a, coe, P, x) {

	if (is.vector(x)) {
        	x <- matrix(x)
    	}
	if (is.vector(P)) {
        	P <- matrix(rep(P,each=ncol(x)),ncol=nrow(x),byrow=TRUE)
    	}

	mu<-a


	if(mu<0){
		mu=1e-22
	}
	if(mu>1){
		mu=1-1e-22
	}

	mu<-rep(mu,ncol(P))

	c<-P[1,]
	d<-rep(1,length(c))
	p<-P[2,]
	o<-P[3,]

    	f <- 0

	xi<-as.numeric(x[,3])
	ni<-as.numeric(x[,4])+xi

	s<-x[,2]

	f<-f+funcv.c("i",s, mu, xi, ni, c,d,p,o)

    return(-f)
	
}

.opl.u.i.c.d <- function(a, P, x) {

	if (is.vector(x)) {
        	x <- matrix(x)
    	}
	if (is.vector(P)) {
        	P <- matrix(rep(P,each=ncol(x)),ncol=nrow(x),byrow=TRUE)
    	}

	c<-rep(a,ncol(P))

	mu<-P[1,]
	d<-rep(1,length(c))
	p<-P[2,]
	o<-P[3,]

    	f <- 0

	xi<-as.numeric(x[,3])
	ni<-as.numeric(x[,4])+xi

	s<-x[,2]

	f<-f+funcv.c("i",s, mu, xi, ni, c,d,p,o)

    return(-f)
	
}

.opl.u.c.c <- function(a, P, x, index) {

	if (is.vector(x)) {
        	x <- matrix(x)
    	}
	if (is.vector(P)) {
        	P <- matrix(rep(P,each=ncol(x)),ncol=ncol(x),byrow=TRUE)
    	}

	c<-a[index]

	mu<-P[1,]

	d<-rep(1,length(c))
	p<-P[2,]
	o<-P[3,]

    	f <- 0

	xi<-as.numeric(x[,3])
	ni<-as.numeric(x[,4])+xi

	s<-x[,2]

	f<-f+funcv.c("du",s, mu, xi, ni, c,d,p,o)

    return(-f)
	
}
.opl.g.c.c <- function(a, P, x, index) {

	if (is.vector(x)) {
        	x <- matrix(x)
    	}
	if (is.vector(P)) {
        	P <- matrix(rep(P,each=ncol(x)),ncol=ncol(x),byrow=TRUE)
    	}

	c<-a[index]

	mu<-P[1,]

	d<-rep(1,length(c))
	p<-P[2,]
	o<-P[3,]

    	g <- rep(0,length(a))

	xi<-as.numeric(x[,3])
	ni<-as.numeric(x[,4])+xi

	s<-x[,2]

	g<-g+funcv.c("dg",s, mu, xi, ni, c,d,p,o,index)

    return(-g)
	
}


.opl.g.coe.fix.c <- function(P, coe, group, pi.m, ref, len.r, x) {

    	if (is.vector(x)) {
        	x <- matrix(x)
    	}

	c<-P[1]
	d<-1
	p<-P[2]
	o<-P[3]

	if(c<=0){
		c=1e-5
	}
	if(p<0){
		p=0
	}
	if(p>1){
		p=1
	}
	if(o<0){
		o=0
	}

    	g <- rep(0,length(P))

	for (m in 1:(length(group)-1)){
		if(group[m]=="N"){
			next
		}
		for(n in (m+1):length(group)){
			if(group[n]=="N"){
				next
			}
			if(group[m]==group[n]){
				next
			}

			mu<-list()
			for(r in ref){
				mu[[r]] <- rep(pi.m[[r]][m,n],len.r[[r]])
			}
			mu<-unlist(mu)

			xi<-as.numeric(x[,m+2])
			ni<-as.numeric(x[,n+2])+xi

			s<-as.numeric(x[,2])

			g<-g+funcv.c("g", s, mu, xi, ni, c,d,p,o)

			
		}							

	}

    	return(-g)
}

.opl.u.coe.fix.c <- function(P, coe, group, pi.m, ref, len.r, x) {

	if (is.vector(x)) {
        	x <- matrix(x)
    	}

	c<-P[1]
	d<-1
	p<-P[2]
	o<-P[3]

	if(c<=0){
		c=1e-5
	}
	if(p<0){
		p=0
	}
	if(p>1){
		p=1
	}
	if(o<0){
		o=0
	}

    	f <- 0

	for (m in 1:(length(group)-1)){
		if(group[m]=="N"){
			next
		}
		for(n in (m+1):length(group)){
			if(group[n]=="N"){
				next
			}
			if(group[m]==group[n]){
				next
			}

			mu<-list()
			for(r in ref){
				mu[[r]] <- rep(pi.m[[r]][m,n],len.r[[r]])
			}
			mu<-as.numeric(unlist(mu))

			xi<-as.numeric(x[,m+2])
			ni<-as.numeric(x[,n+2])+xi

			s<-as.numeric(x[,2])

			f<-f+funcv.c("u",s, mu, xi, ni, c,d,p,o)

			
		}

	}

	return(-f)
}


.opl.u.i <- function(a, coe, P, x) {

	if (is.vector(x)) {
        	x <- matrix(x)
    	}

	mu<-a


	if(mu<0){
		mu=1e-22
	}
	if(mu>1){
		mu=1-1e-22
	}

	c<-P[1]
	d<-1
	p<-P[2]
	o<-P[3]

    	f <- 0

	for (i in 1:nrow(x)) {

		ref<-x[i,1]
		strand<-x[i,2]

		coe.s<-coe[[strand]]

    		word.m<-unlist(x[i,3:nrow(coe.s)])
		word<-which(word.m==1)
		D<-c*exp(d*sum(coe.s[word,2]))

		count<-as.numeric(unlist(x[i,(nrow(coe.s)+3):(ncol(x))]))

				xi <- count[1]
        			ni <- count[2]+xi
		
				theta<-D/((ni+o)^p)

        			if (xi>0) {
            				f <- f + log(mu)

            				if (xi>1) {
                				r <- 1:(xi-1)
                				v <- mu/r+theta
                				f <- f + sum(log(v))
            				}
        			}

        			if (xi < ni) {
            				f <- f + log(1-mu)
            				if (xi < (ni-1)) {
                				r <- 1:(ni-xi-1)
                				v <- (1-mu)/r+theta
                				f <- f + sum(log(v))
            				}
        			}

        			if (ni > 1) {
            				r <- 1:(ni-1)
            				v <- 1/r + theta
            				f <- f - sum(log(v))
        			}
	}

	return(-f)
}


.opl.g.coe.fix <- function(P, coe, group, pi.m, x) {

    	if (is.vector(x)) {
        	x <- matrix(x)
    	}

	c<-P[1]
	d<-1
	p<-P[2]
	o<-P[3]

	if(c<=0){
		c=1e-5
	}
	if(p<0){
		p=0
	}
	if(p>1){
		p=1
	}
	if(o<0){
		o=0
	}

    	g <- rep(0,length(P))

    	for (i in 1:nrow(x)) {

		ref<-x[i,1]
		strand<-x[i,2]

		coe.s<-coe[[strand]]

    		word.m<-unlist(x[i,3:(nrow(coe.s)+2)])
		word<-which(word.m==1)
		D<-c*exp(d*sum(coe.s[word,2]))

		count<-as.numeric(unlist(x[i,(nrow(coe.s)+3):(ncol(x))]))

		for (m in 1:(length(group)-1)){
			if(group[m]=="N"){
				next
			}
			for(n in (m+1):length(group)){
				if(group[n]=="N"){
					next
				}
				if(group[m]==group[n]){
					next
				}

				mu <- pi.m[[ref]][m,n]

				xi <- count[m]
        			ni <- count[n]+xi

				if(ni<1){
					next;
				}
		
				theta<-D/((ni+o)^p)

				gt<-0

				if (xi>0) {

            				if (xi>1) {
                				r <- 1:(xi-1)
                				v <- mu/r+theta

                				one_on_v <- 1 / v

						gt <- gt + sum(one_on_v)
            				}
        			}

        			if (xi < ni) {

            				if (xi < (ni-1)) {
                				r <- 1:(ni-xi-1)
                				v <- (1-mu)/r+theta
                				one_on_v <- 1 / v

						gt <- gt + sum(one_on_v)

            				}
        			}

        			if (ni > 1) {
            				r <- 1:(ni-1)
            				v <- 1/r + theta
             				one_on_v <- 1 / v
            				gt <- gt - sum(one_on_v)
        			}

				g[1]<-g[1]+gt*theta/c
				g[2]=g[2]+gt*(-1)*D*log(ni+o)*(ni+o)^(-p)
				g[3]=g[3]+gt*(-1)*p*D*(ni+o)^(-p-1)

			}
		} 							

	}

    	return(-g)
}

.opl.u.coe.fix <- function(P, coe, group, pi.m, x) {

	if (is.vector(x)) {
        	x <- matrix(x)
    	}

	c<-P[1]
	d<-1
	p<-P[2]
	o<-P[3]

	if(c<=0){
		c=1e-5
	}
	if(p<0){
		p=0
	}
	if(p>1){
		p=1
	}
	if(o<0){
		o=0
	}

    	f <- 0

	for (i in 1:nrow(x)) {

		ref<-x[i,1]
		strand<-x[i,2]

		coe.s<-coe[[strand]]

    		word.m<-unlist(x[i,3:(nrow(coe.s)+2)])
		word<-which(word.m==1)
		D<-c*exp(d*sum(coe.s[word,2]))

		count<-as.numeric(unlist(x[i,(nrow(coe.s)+3):(ncol(x))]))

		for (m in 1:(length(group)-1)){
			if(group[m]=="N"){
				next
			}
			for(n in (m+1):length(group)){
				if(group[n]=="N"){
					next
				}
				if(group[m]==group[n]){
					next
				}

				mu <- pi.m[[ref]][m,n]

				xi <- count[m]
        			ni <- count[n]+xi

				if(ni<1){
					next
				}
		
				theta<-D/((ni+o)^p)

        			if (xi>0) {
            				f <- f + log(mu)

            				if (xi>1) {
                				r <- 1:(xi-1)
                				v <- mu/r+theta
                				f <- f + sum(log(v))
            				}
        			}

        			if (xi < ni) {
            				f <- f + log(1-mu)
            				if (xi < (ni-1)) {
                				r <- 1:(ni-xi-1)
                				v <- (1-mu)/r+theta
                				f <- f + sum(log(v))
            				}
        			}

        			if (ni > 1) {
            				r <- 1:(ni-1)
            				v <- 1/r + theta
            				f <- f - sum(log(v))
        			}
			}
		}
	}

	return(-f)
}

.opl.u.f.c <- function(a, p, o, x, u) {

	if (is.vector(x)) {
        	x <- matrix(x)
    	}

	p<-p
	o<-o

    	f <- 0

	mu<-u

    	for (i in 1:nrow(x)) {

		word<-unlist(x[i,1:13])
		D<-exp(sum(a[word]))

        	xi <- x[i,14]
        	ni <- x[i,15]+xi
		

		theta<-D/((ni+o)^p)


        	if (xi>0) {
            		f <- f + log(mu)

            		if (xi>1) {
                		r <- 1:(xi-1)
                		v <- mu/r+theta
                		f <- f + sum(log(v))
            		}
        	}

        	if (xi < ni) {
            		f <- f + log(1-mu)
            		if (xi < (ni-1)) {
                		r <- 1:(ni-xi-1)
                		v <- (1-mu)/r+theta
                		f <- f + sum(log(v))
            		}
        	}

        	if (ni > 1) {
            		r <- 1:(ni-1)
            		v <- 1/r + theta
            		f <- f - sum(log(v))
        	}

    	}

    	return(-f)
}

.opl.u.c.One <- function(a, o, x) {

	if (is.vector(x)) {
        	x <- matrix(x)
    	}

	p<-a[length(a)]
	o<-o

    	f <- 0

	for(i in 1:nrow(x)){

		mu <- x[i,ncol(x)]

    		word<-unlist(x[i,1:((length(a)-1)/4)])
		D<-exp(sum(a[word]))

        	xi <- x[i,ncol(x)-2]
        	ni <- x[i,ncol(x)-1]+xi
		

		theta<-D/((ni+o)^p)

        	if (xi>0) {
            		f <- f + log(mu)

            		if (xi>1) {
                		r <- 1:(xi-1)
                		v <- mu/r+theta
                		f <- f + sum(log(v))
            		}
        	}

        	if (xi < ni) {
            		f <- f + log(1-mu)
            		if (xi < (ni-1)) {
                		r <- 1:(ni-xi-1)
                		v <- (1-mu)/r+theta
                		f <- f + sum(log(v))
            		}
        	}

        	if (ni > 1) {
            		r <- 1:(ni-1)
            		v <- 1/r + theta
            		f <- f - sum(log(v))
        	}
	}

	return(-f)
}

.opl.g.c.D <- function(a, x) {

    	if (is.vector(x)) {
        	x <- matrix(x)
    	}

	p<-a[length(a)-1]
	o <-a[length(a)]

    	g <- rep(0,length(a))

    	for (i in 1:nrow(x)) {

        	mu <- x[i,ncol(x)]

    		word<-unlist(x[i,1:((length(a)-2)/4)])
		D<-exp(sum(a[word]))

        	xi <- x[i,ncol(x)-2]
        	ni <- x[i,ncol(x)-1]+xi
		

		theta<-D/((ni+o)^p)

		gt<-0

		if (xi>0) {

            		if (xi>1) {
                		r <- 1:(xi-1)
                		v <- mu/r+theta

                		one_on_v <- 1 / v

				gt <- gt + sum(one_on_v)
            		}
        	}

        	if (xi < ni) {

            		if (xi < (ni-1)) {
                		r <- 1:(ni-xi-1)
                		v <- (1-mu)/r+theta
                		one_on_v <- 1 / v

				gt <- gt + sum(one_on_v)

            		}
        	}

        	if (ni > 1) {
            		r <- 1:(ni-1)
            		v <- 1/r + theta
             		one_on_v <- 1 / v
            		gt <- gt - sum(one_on_v)
        	}

		g[length(a)-1]=g[length(a)-1]+gt*(-1)*D*log(ni+o)*(ni+o)^(-p)
		g[length(a)]=g[length(a)]+gt*(-1)*p*D*(ni+o)^(-p-1)

		gt<-gt*theta
			
		g[word]<-g[word]+gt			

	}

    	return(-g)
}

.opl.u.c.D <- function(a, x) {

	if (is.vector(x)) {
        	x <- matrix(x)
    	}

	p<-a[length(a)-1]
	o<-a[length(a)]

    	f <- 0

	value<-list()

	for(i in 1:nrow(x)){

		mu <- x[i,ncol(x)]

    		word<-unlist(x[i,1:((length(a)-2)/4)])
		D<-exp(sum(a[word]))
		

        	xi <- x[i,ncol(x)-2]
        	ni <- x[i,ncol(x)-1]+xi
		

		theta<-D/((ni+o)^p)

        	if (xi>0) {
            		f <- f + log(mu)

            		if (xi>1) {
                		r <- 1:(xi-1)
                		v <- mu/r+theta
                		f <- f + sum(log(v))
            		}
        	}

        	if (xi < ni) {
            		f <- f + log(1-mu)
            		if (xi < (ni-1)) {
                		r <- 1:(ni-xi-1)
                		v <- (1-mu)/r+theta
                		f <- f + sum(log(v))
            		}
        	}

        	if (ni > 1) {
            		r <- 1:(ni-1)
            		v <- 1/r + theta
            		f <- f - sum(log(v))
        	}
		value[[i]]<-f
	}

	return(-f)
}


x.test <- function(x,pi.i,pi.bi.i,pn.i,a0) {
	
	fval.bb.pi<-.opl.u.i.c(pi.i,coe=NULL,a0,x)*(-1)
	fval.bb.pn<-.opl.u.i.c(pn.i,coe=NULL,a0,x)*(-1)
	p.bb<-.bb.chis.test(fval.bb.pi,fval.bb.pn)

	fval.bi.pi<-.bi.f(rep(pi.bi.i,nrow(x)),x[,3:4])*(-1)
	fval.bi.pn<-.bi.f(rep(pn.i,nrow(x)),x[,3:4])*(-1)
	p.bi<-.bb.chis.test(fval.bi.pi,fval.bi.pn)

	return(c(p.bb,p.bi))
}

.bb.chis.test <- function(fval.pi,fval.pn,df=1) {

	B <- 2*(fval.pi - fval.pn)
        p = pchisq(B, df = df, lower.tail=FALSE)
	
	if(p<1e-300){
		p<-0
	}

	return(p)
}

.bi.f <- function(a,x) {

	if (is.vector(x)) {
        	x <- matrix(x, nrow=1)
    	}

	pi.ori<-a

	bi.val<-0
	for(i in 1:length(x[,1])){
		xi<-as.numeric(x[i,1])
		ni<-as.numeric(x[i,2])+xi

		pi<-pi.ori[i]

		bi.val<-bi.val+dbinom(xi,ni,pi,log=TRUE)
	}

        return(-bi.val)
}
