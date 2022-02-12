MY_SSH=$1
dx download -r $MY_SSH
chmod 400 .ssh/id_ed25519
eval "$(ssh-agent -s)"
ssh-add .ssh/id_ed25519
git clone git@github.com:jsture/gogoGPCR.git
cd gogoGPCR/
pip install -e .

