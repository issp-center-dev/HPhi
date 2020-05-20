# This is a pull request, finish.
if [ "_$TRAVIS_PULL_REQUEST" != "_false" ] ;then
  echo "This is a pull request, do nothing."
  exit 0;
fi
# build doc if and only if master, develop, xxx-autodoc, and tag
feature_branch=${TRAVIS_BRANCH%-autodoc}

if [ "_$TRAVIS_BRANCH" == "_master" ]; then
  echo "This is the master branch, deploy docs."
elif [ "_$TRAVIS_BRANCH" == "_develop" ]; then
  echo "This is the develop branch, deploy docs."
elif [ "_${feature_branch}" != "_${TRAVIS_BRANCH}" ]; then
  echo "This is an auto-documented branch, deploy docs."
elif [ -n "$TRAVIS_TAG" ]; then
  echo "This is a versioned tag, deploy docs."
else
  echo "Do nothing."
  exit 0
fi

set -e

sudo apt-get install -y enchant
sudo apt-get install -y texlive-latex-recommended texlive-latex-extra texlive-lang-japanese texlive-fonts-recommended texlive-fonts-extra latexmk
kanji-config-updmap-sys ipaex
sudo apt-get install -y python3 python3-sphinx python3-sphinxcontrib.spelling

openssl aes-256-cbc -K $encrypted_87f43018402c_key -iv $encrypted_87f43018402c_iv -in ${ROOTDIR}/.travis_scripts/id_rsa.enc -out ~/.ssh/id_rsa -d

chmod 600 ~/.ssh/id_rsa
echo -e "Host github.com\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config

git clone git@github.com:${TRAVIS_REPO_SLUG} hphi-doc
cd hphi-doc
mkdir build && cd build
cmake -DDocument=ON ../
make doc-ja-html
make doc-en-html

set +e

git checkout gh-pages
if [ ${feature_branch} != ${TRAVIS_BRANCH} ]; then
  docdir=${feature_branch}
elif [ "_${TRAVIS_BRANCH}" == "_develop" ]; then
  docdir=develop
elif [ "_${TRAVIS_BRANCH}" == "_master" ]; then
  docdir=master
elif [ -n ${TRAVIS_TAG}]; then
  docdir=${TRAVIS_TAG}
else
  echo "The deploy script failed to solve where to install documents. The script has some mistake."
  echo "\$TRAVIS_BRANCH: $TRAVIS_BRANCH"
  echo "\$TRAVIS_TAG: $TRAVIS_TAG"
  echo "\$TRAVIS_PULL_REQUEST: $TRAVIS_PULL_REQUEST"
  echo "\$feature_branch: $feature_branch"
  exit 1
fi

cd ${ROOTDIR}/hphi-doc/manual
mkdir -p $docdir && cd $docdir
for lang in ja en; do
  rm -rf $lang/html
  mkdir -p $lang
  cp -r ${ROOTDIR}/hphi-doc/build/doc/${lang}/source/html $lang/html
  git add $lang/html
done

git config --global user.email "hphi-dev@issp.u-tokyo.ac.jp"
git config --global user.name "HPhi"
git commit -m "Update by TravisCI"
ST=$?
if [ $ST == 0 ]; then
  git push origin gh-pages:gh-pages --follow-tags > /dev/null 2>&1
fi

