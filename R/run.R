scp -r www/ rathik@reslnmarisws01.research.chop.edu:"/data/transport/IMMUNO_PROFILE/"
ssh rathik@reslnmarisws01.research.chop.edu
rsync -avhH /data/transport/IMMUNO_PROFILE rathik@reslnmaris01:/home/rathik
ssh rathik@reslnmaris01
rsync -avhH /home/rathik/IMMUNO_PROFILE /data/shiny-server