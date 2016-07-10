#include "mex.h"
#include <math.h>
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]){
    
    int *cmax,  imaxPre, tots, imax, curri, t, i, k,  inner_it, nt0, NT, type, maxiter, *inner_all;
    double *c, *Params, *st, *kernel, *F0, *WtW, *Wf, *Mmax, Th, Thi, M0max;
    double  *delta, ncoef, *new_coef, newD, *DD;
        
    NT      = (int) mxGetM(prhs[1]);
    nt0     = (int) mxGetM(prhs[2]);
    
   //mexPrintf("i exist\n");

    Params  = mxGetPr(prhs[0]); /* 1 by 1 */
    type    = (int) Params[0];
    Th      = (double) Params[1];
    Thi     = (double) Params[2];
    maxiter = (int) Params[3];
    F0      = mxGetPr(prhs[1]); /* 1 by 1 */
    kernel  = mxGetPr(prhs[2]); /* 1 by 1 */
    
    
    plhs[0]  = mxCreateDoubleMatrix(maxiter, 1, mxREAL);
    c        = mxGetPr(plhs[0]);
    plhs[1]  = mxCreateDoubleMatrix(maxiter, 1, mxREAL);
    st        = mxGetPr(plhs[1]);
    
//    plhs[1]  = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
 //   st       = (int*) mxGetData(plhs[1]);
//    plhs[2]  = mxCreateNumericArray(1, maxiter, mxINT32_CLASS, mxREAL);
 //   inner_all= (int*) mxGetData(plhs[2]);
    
        
    WtW      = (double*) calloc(2*nt0-1, sizeof(double));
    Wf       = (double*) calloc(NT,      sizeof(double));
    delta    = (double*) calloc(maxiter, sizeof(double));
    Mmax     = (double*) calloc(maxiter, sizeof(double));
    cmax     = (int*)    calloc(maxiter, sizeof(int));
    DD       = (double*) calloc(maxiter*(2*nt0-1), sizeof(double));
    new_coef = (double*) calloc(maxiter*(2*nt0-1), sizeof(double));
    

    //compute WtW
    for (t=nt0;t<2*nt0;t++)
        for (k=0;k<nt0;k++)
            if (t-k-1>=0 && t-k-1<2*nt0-1)
                WtW[t-k-1] = WtW[t-k-1] + kernel[t-nt0] * kernel[k];
    
    //compute Wf 
    for (t=0;t<NT;t++)
        for (k=0;k<nt0;k++)
            if (t-k>=0 && t-k<NT)
                Wf[t-k] = Wf[t-k] + F0[t] * kernel[k];
    
    
    //big loop   
    tots = 0;
    while (tots<maxiter){
        inner_it = 0;      
     //   mexPrintf("iteration %d\n", tots);
         // small loop        
        while (tots>0) {
            M0max = 0.0f;
            imax  = 0;
            // compute stats for best new position
            for (i=0;i<tots;i++){
                if ((inner_it==0 && abs(st[tots-1]-st[i])<2*nt0) ||
                        (inner_it>0 && abs(st[imaxPre]-st[i])<2*nt0)){
                    // estimate delta loss upon discarding this spike
                    delta[i] = c[i] * c[i] + 2 * c[i] * Wf[(int) st[i]];
                    Mmax[i]  = 0;
                    cmax[i]  = 0;
                    for (t=0; t<2*nt0-1;t++){
                        curri = st[i] + t+1 - nt0;
                        if (curri>=0 && curri<NT){
                            ncoef                     = Wf[curri] + WtW[t] * c[i];
                            new_coef[t + i*(2*nt0-1)] = ncoef;
                            // estimate increase in explained variance from shifting each spike
                            newD                      = max(0, ncoef) * max(0, ncoef) - delta[i];
                            DD[t + i * (2*nt0-1)]     = newD;
                            // find best shift for each spike
                            if (newD > Mmax[i]){
                                Mmax[i] = newD;
                                cmax[i] = t;
                            }
                        }
                    }
                    
                }
                // find best shift overall
                if (Mmax[i]>M0max){
                    M0max = Mmax[i];
                    imax  = i;
                }
                //}
            }
            
//             if (inner_it==0){
//                 c[0] = M0max;
//                 //c[0] = Wf[(int) st[i]];
//                 //c[0] = WtW[nt0-1];
//                 // c[0] = cmax[0];
//                 break;
//             }
            
            if (M0max<Thi){
                // exit the inner loop if it's not worth it
                break;
            }
            else{
                // add back contribution of shifted spike
                for (t=0; t<2*nt0-1;t++){
                    curri      = st[imax] + t -(nt0-1);
                    if (curri>=0 && curri<NT)
                        Wf[curri] += c[imax] * WtW[t];
                }
                for (t=0; t<nt0;t++){
                    curri = st[imax] + t;
                    if (curri<NT)
                        F0[curri]        += c[imax] * kernel[t]; 
                }
                
                // re-assign spike time and magnitude
                c[imax]  = new_coef[cmax[imax] + imax * (2*nt0-1)];
                st[imax] += cmax[imax] - (nt0-1);      
                imaxPre = imax;
                
                // remove contribution of shifted spike
               for (t=0; t<2*nt0-1;t++){
                    curri      = st[imax] + t -(nt0-1);
                    if (curri>=0 && curri<NT)
                        Wf[curri] -= c[imax] * WtW[t];
                }
                for (t=0; t<nt0;t++){
                    curri = st[imax] + t;
                    if (curri<NT)
                        F0[curri]        -= c[imax] * kernel[t]; 
                }
             
                 
                // increment counter of spike reshufflings
                inner_it++;
            }
                
        }
         
       // inner_all[tots] = inner_it;

        
        // find the new best place to introduce a spike
        M0max = 0.0f;
        for (t=0;t<NT;t++)
            if (Wf[t]>M0max){
                M0max = Wf[t];
                imax = t;
            }
        if (M0max<Th) break;
                
        
        c[tots]   = M0max;
        st[tots]  = imax;
        
        // subtract kernel from raw data
        for(i=0;i<nt0;i++){
            curri = imax + i;
            if (curri<NT)
                F0[imax + i] = F0[imax + i] - M0max * kernel[i];
        }
        
        // subtract kernel from filtered values
        for(i=0;i<2*nt0-1;i++){
            curri = imax + i-(nt0-1);
            if (curri>=0 && curri<NT)
                Wf[curri] = Wf[curri] - M0max * WtW[i];
        }
        tots++;
    }
    
    free(WtW); free(Wf); free(delta); free(cmax); free(DD); free(Mmax); free(new_coef);
    return;
}