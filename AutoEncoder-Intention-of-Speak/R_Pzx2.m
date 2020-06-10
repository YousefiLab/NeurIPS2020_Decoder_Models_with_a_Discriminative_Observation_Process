function [Pzx,pk,qk,zk]=R_Pzx2(Data,Par,xvals,ii)
% compute p(zk|xk)
ipk=find(Data.train.ti(ii) - Data.ZRes< Data.train.pwords & Data.train.pwords <= Data.train.ti(ii));
iqk=find(Data.train.ti(ii) - Data.ZRes< Data.train.cwords & Data.train.cwords <= Data.train.ti(ii));

ipkc=find(Data.train.ti(ii-1) - Data.ZHis< Data.train.pwords & Data.train.pwords <= Data.train.ti(ii-1));
iqkc=find(Data.train.ti(ii-1) - Data.ZHis< Data.train.cwords & Data.train.cwords <= Data.train.ti(ii-1));
        if ~isempty(ipk) && ~isempty(iqk)
            zk=1;
            pk= length(ipkc);
            qk= length(iqkc);
           Landa=exp(Par.Z.a0 + Par.Z.a1 * xvals + Par.Z.a2*pk + Par.Z.a3*qk);
           Pzx= Landa;
           if Pzx ~=0
           Pzx=Pzx./sum(Pzx); 
           else
           Pzx=1./length(Pzx);     
           end
        elseif  isempty(ipk) && ~isempty(iqk)
            zk=0;
            pk= 0;
            qk= length(iqkc);
            Landa=exp(Par.Z.a0 + Par.Z.a1 * xvals  + Par.Z.a3*qk);
           Pzx=ones(size(Landa));
           if Pzx ~=0
           Pzx=Pzx./sum(Pzx); 
           
           end
           
        elseif  ~isempty(ipk) && isempty(iqk)
            zk=1;
            pk= length(ipkc);
            qk= 0;
            Landa=exp(Par.Z.a0 + Par.Z.a1 * xvals  + Par.Z.a2*pk);
           Pzx= Landa;
           
           if Pzx ~=0
           Pzx=Pzx./sum(Pzx); 
           
           end
        else
            pk=0;
            qk= 0;
            zk=0;
           Landa=exp(Par.Z.a0 + Par.Z.a1 * xvals  );
           Pzx=ones(size(Landa));
           if Pzx ~=0
           Pzx=Pzx./sum(Pzx); 
           
           end
        end
