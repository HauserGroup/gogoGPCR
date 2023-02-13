MY_SSH=$1
GIT_EMAIL=$2
GIT_USER=$3
dx download -r $MY_SSH
chmod 400 .ssh/id_ed25519
eval "$(ssh-agent -s)"
ssh-add .ssh/id_ed25519
git config --global user.email $GIT_EMAIL
git config --global user.name $GIT_USER
git clone git@github.com:HauserGroup/gogoGPCR.git
cd gogoGPCR/
pip install -e .

