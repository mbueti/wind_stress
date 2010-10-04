% Sample MATLAB code to read and plot wind stress

im=508;
jm=449;

latmin=10.;
latmax=47.5;
longmin=-98.5 ;
longmax=-50.;

var='TXY';
dates=[{'05082600'},{'05082606'},{'05082612'},{'05082618'},{'05082700'},{'05082706'},{'05082712'},{'05082718'},{'05082800'},{'05082806'},{'05082812'},{'05082818'},{'05082900'},{'05082906'},{'05082912'},{'05082918'},{'05083000'},{'05083006'},{'05083012'}];

x=zeros(im);
y=zeros(jm);
for j=1:jm
  for i=1:im
    y(j)=latmin+(j-1)*(latmax-latmin)/(jm-1);
    x(i)=longmin+(i-1)*(longmax-longmin)/(im-1);
  end
end

TX=zeros(im,jm,size(dates,2));
TY=zeros(im,jm,size(dates,2));
for i=1:size(dates,2)
  date=cell2mat(dates(i));
  filename=['../output/',var,'.',date];
  fid=fopen(filename,'r','b');
  n=fread(fid,1,'int32');
  skip=0;
  fseek(fid,skip,'cof');
  TX(:,:,i)=fread(fid,[im,jm],'float');
  n=fread(fid,1,'int32');
  n=fread(fid,1,'int32');
  TY(:,:,i)=fread(fid,[im,jm],'float');
  fclose(fid);
end


TXY=sqrt(TX.^2+TY.^2);

%pcolor(x,y,TXY')
%shading interp
%xlim([-90 -85])
%ylim([23 28])
%print -dpng windstress
%close all
