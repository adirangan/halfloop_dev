rsync -avum \
    /data/rangan/dir_bcc/dir_halfloop_dev \
    /home/rangan/dir_bcc/ \
    --include="*/" \
    --include="*.h" \
    --include="*.c" \
    --include="*.in" \
    --include="*.make" \
    --include="*.sh" \
    --exclude="*" ;
cd /data/rangan/dir_bcc/dir_halfloop_dev/ ;
git add dir_h/*.h ;
git add dir_c/*.c ;
git add dir_in/*.in ;
git add *.make ;
git add *.sh ;
git commit -m "updating dir_halfloop_dev with halfloop_split " ;
git push ;
git pull ;
cd /data/rangan/dir_bcc/dir_halfloop_dev ;

#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
#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
#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
#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 
#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
#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
#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
#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
#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
#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
#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
#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
