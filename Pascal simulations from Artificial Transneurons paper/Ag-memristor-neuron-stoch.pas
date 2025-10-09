Program simul;
uses math;
 var temp, dt,  tmax, t, r, i0,i00, i02, imax, avervolt, ttt0,
 kap, lam,x0,avr,avoverr,avtemp, Tb, chx, output, output1,output2, xx, xxxx,ppot,fforce, 
 
 r0,tauc, rmin, ccc:real;

        lll,l,mmm,aver,mmt:longint;

        x: array[0..201] of real;
        vx: array[1..200] of real;
             
        res,res1, resx,respot: text;

        label  ord1;

         function fT(time:extended):real;
         var rana,R1,R2,XR,YR:real;
         begin
         rana:=random;
         R1:=-ln(1.0-rana);
         rana:=random;
         R2:=2*pi*rana;
         R1:=sqrt(2.0*R1);
         XR:=R1*cos(R2);
         YR:=R1*sin(R2);
         fT:=sqrt(2*temp/dt)*XR;
         end;
        
         
         function pot(xxx:real):real;
                
         
         begin
                                
         pot:=-(sqr(xxx-0.1)+0.1)-180*exp(-sqr(xxx-0.8)/(0.01))+0.2*sqrt(10)*i0
                  
         end;
        
        
            begin
            aver:=10000;
       
            assign(res1,'info-bet-2d.dat'); assign(respot,'pot-bet-2d.dat');             
            rewrite(res1); rewrite(respot);             
            assign(resx,'x.dat');
            rewrite(resx);
            assign(res,'ri-bet-2d1.dat'); 
            rewrite(res); 

            writeln('number of Ag-clusters'); readln(mmm); mmt:=1; writeln(res1, 'mmm=', mmm);
            
            x[0]:=-1; x[mmm+1]:=1;
            
            randomize;
             
           
            writeln('kappa'); readln(kap); writeln(res1,'kappa=',kap);
            x0:=0.85; writeln(res1,'x0=',x0);
            writeln('Tb'); readln(Tb); writeln(res1,'Tb=',Tb);
            writeln('maximum voltage imax'); readln(imax); writeln(res1,'imax',imax);
            tauc:=30; writeln(res1,'tauc=',tauc); 
            writeln('r0'); readln(r0); writeln(res1,'r0=',r0);
            lam:=0.12; writeln(res1, 'lam=', lam);
                       
            i0:=imax; writeln('ccc'); readln(ccc); writeln(res1,'ccc=',ccc);
            xxxx:=-1.1; ppot:=0; fforce:=0;  while xxxx<=1.1 do begin fforce:=pot(xxxx); ppot:= ppot-fforce*0.0001;
            
            writeln(respot,xxxx,';',fforce,';',ppot); xxxx:=xxxx+0.0001 end; 
            i0:=0; 
           writeln(respot); writeln(respot);       
           xxxx:=-1.1; ppot:=0; fforce:=0;  while xxxx<=1.1 do begin fforce:=pot(xxxx); ppot:= ppot-fforce*0.0001; writeln(respot,xxxx,';',fforce,';',ppot); 
           xxxx:=xxxx+0.0001 end; readln;
           
            close(respot); 
           
            l:=1; while l<=mmm do begin
            x[l]:=-1+0.001*l; 
            l:=l+1 end;
                      
            i0:=0; 
            
            temp:=Tb; ttt0:=Tb; t:=0;
           
                   
            tmax:=50000;
            dt:=0.00001; 
                    
            close(res1);
            
            
            lll:=1; avr:=0; avoverr:=0; avtemp:=0; avervolt:=0; 
                                    
                                        
            while t<tmax do begin
                                    
            i00:=imax; 
             
                      
            if mmm = 2 then begin l:=1; r:=0; while l<=mmm+1 do begin
            r:=r+exp((x[l]-x[l-1])/lam); l:=l+1 end; r:=r/((mmt+1)*exp(2/((mmt+1)*lam))); end;
           
            if mmm=1 then begin r:=0.5*(exp(x[1]/lam)+exp(-x[1]/lam)) end;

            i02:=i0*i0; 
            
            temp:=temp+(ccc*i02/r-kap*(temp-ttt0))*dt;
                                
            if mmm = 2 then begin vx[2]:=ft(t)+pot(x[2]); vx[1]:=ft(t)+pot(x[1]) end; 
            if mmm=1 then vx[1]:=ft(t)+pot(x[1]);
            
            l:=1; while l<=mmm do begin
            x[l]:=x[l]+vx[l]*dt; l:=l+1;
            end;

            i0:=i0+(dt/tauc)*(i00-((r0/r)+1)*i0);  
            
            l:=1; while l<=mmm do begin
            if abs(x[l])>1 then begin x[l]:=x[l]/abs(x[l])-0.000001*x[l]/abs(x[l]) end; l:=l+1 end;
          
            ord1: l:=1; while l<=mmm-1 do begin if x[l]>x[l+1] then begin chx:=x[l]; x[l]:=x[l+1]; x[l+1]:=chx; goto ord1 end; l:=l+1 end;
                              
            avr:=avr+r; avoverr:=avoverr+1/r; avtemp:=avtemp+temp;           
            avervolt:=avervolt+i0; 
               
                        
            if frac(lll/100000)=0 then begin writeln(' t=',t) end;
            if frac(lll/aver)=0 then begin output:=avr/aver; output1:=avoverr/aver; output2:=avtemp/aver; avervolt:=avervolt/aver; 
            writeln(res,t,';',i00,';',i0,';',output,';',output1,';',output2,';',ttt0,';',avervolt);
           
            write(resx,t,' ',x[0]); l:=1; while l<=mmm+1 do begin
            write(resx,' ',x[l]); l:=l+1 end; writeln(resx);

            avr:=0; avoverr:=0; avtemp:=0; avervolt:=0; end; t:=t+dt; lll:=lll+1;

           end;    
          
 close(res); close(resx); 

            end.

