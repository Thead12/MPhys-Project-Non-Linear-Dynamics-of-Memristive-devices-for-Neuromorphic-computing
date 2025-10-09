Program simul;
uses math;
 var dt,  tmax, t, r, i0,i00, imax, avervolt, 
 kap, lam,x0,avr,avoverr, output, output1, xxxx,ppot,fforce, tempforce, tempforce1, dr, qtemp,
 
 r0,tauc, ccc:real;

        lll,l,mmm,aver:longint;

        x: array[0..201] of real;
        vx: array[1..200] of real;
             
        res,res1, resx,respot: text;

      
      
         
         function pot(xxx:real):real;
                
         
         begin
                                
         pot:=-(sqr(xxx-0.1)+0.1)-180*exp(-sqr(xxx-0.8)/(0.01))+0.2*sqrt(10)*i0-100*power(xxx,100-1)
                  
         end;
        
        
            begin
            aver:=10000;
            qtemp:=1;
            tempforce:=0;
            tempforce1:=0;
       
            assign(res1,'info-bet-2d.dat'); assign(respot,'pot-bet-2d.dat');             
            rewrite(res1); rewrite(respot);             
            assign(resx,'x.dat');
            rewrite(resx);
            assign(res,'ri-bet-2d1.dat'); 
            rewrite(res); 

            mmm:=1;
            
            x[0]:=-1; x[mmm+1]:=1;
            
            randomize;
             
           
            writeln('kappa'); readln(kap); writeln(res1,'kappa=',kap);
            x0:=0.85; writeln(res1,'x0=',x0);
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
           
            x[1]:=-1+0.001; 
                                  
            i0:=0; 
            
            t:=0;
            qtemp:=1;
            tempforce:=0;
            tempforce1:=0;
                   
            tmax:=50000;
            dt:=0.00001; 
                    
            close(res1);
            
            
            lll:=1; avr:=0; avoverr:=0; avervolt:=0; 
                                    
                                        
            while t<tmax do begin
                                    
            i00:=imax; 
             
                      
            r:=0.5*(exp(x[1]/lam)+exp(-x[1]/lam)); 
            dr:=0.5*(exp(x[1]/lam)-exp(-x[1]/lam))/lam; 

                                  
            vx[1]:=pot(x[1])-qtemp*tempforce;
            tempforce:=tempforce+(ccc*(2*i0*tempforce1*r-i0*i0*dr)/sqr(r)-kap*tempforce)*dt; 
            tempforce1:=tempforce1+(r0*dr*i0/sqr(r)-(r0/r+1)*tempforce1)*dt/tauc;
            i0:=i0+(dt/tauc)*(i00-((r0/r)+1)*i0);  
            
            l:=1; while l<=mmm do begin
            x[l]:=x[l]+vx[l]*dt; l:=l+1;
            end;
                                              
            avr:=avr+r; avoverr:=avoverr+1/r;        
            avervolt:=avervolt+i0; 
                                
            if frac(lll/100000)=0 then begin writeln(' t=',t) end;
            if frac(lll/aver)=0 then begin output:=avr/aver; output1:=avoverr/aver; avervolt:=avervolt/aver; 
            writeln(res,t,';',i00,';',i0,';',output,';',output1,';',avervolt);
           
            write(resx,t,' ',x[0]); l:=1; while l<=mmm+1 do begin
            write(resx,' ',x[l]); l:=l+1 end; writeln(resx);

            avr:=0; avoverr:=0; avervolt:=0; end; t:=t+dt; lll:=lll+1;

           end;    
          
 close(res); close(resx); 

            end.

