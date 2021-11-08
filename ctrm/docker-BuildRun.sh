 #build an image named ctrm
 docker build  -t ctrm  .

#build container and create a shared volume in dp28  at /Shared/CtrmTestData
#u should move input files to this directory
docker run -v /Shared/CtrmTestData:/app/testData  -it ctrm bash 
