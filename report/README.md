Report aimed to be published in BioHackrXiv:
- Style: Markdown converted to PDF with [PDF generator](https://github.com/biohackrxiv/bhxiv-gen-pdf) and around 4 to 10 pages including references
- [Submission guidelines](https://github.com/biohackrxiv/biohackrxiv.github.io/blob/guidelines/submission_guidelines.md)
- [Moderation process](https://github.com/biohackrxiv/biohackrxiv.github.io/blob/guidelines/moderation_process.md)
- [Paper template](https://github.com/biohackrxiv/submission-templates/blob/master/paper.md)

Expected content for mini-publication (taken from mailinglist):
1. Description of problem
2. Description of state-of-the-art
3. Results of work during BioHackathon
4. Future work 

Build paper using docker container (see Docker build instructions [here](https://github.com/biohackrxiv/bhxiv-gen-pdf))
    
    cd $GIT_ROOT
    docker run --rm -it -v $(pwd):/work -w /work biohackrxiv/gen-pdf:local gen-pdf /work/report
