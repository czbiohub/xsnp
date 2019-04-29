#!/bin/bash

set -e

# Usage example:
#
#   devtools/autosave_continuously.sh ubuntu@r:ebs/xsnp
#
# This will monitor for changes and autosave/autopush the current repo both
# to github and to the SSH coordinates above.
#
# Optionally, if there is a Makefile, it will be executed on the SSH host.
#
# When it starts up, it will force-push the current top commit IF THERE
# ARE NO ADDED OR MODIFIED FILES.  This saves any prior work by the user.
#
# Subsequent iterations (one per second) will autosave/autopush/autobuild
# when added or modified files are detected.  The first time an added or
# modified file is detected, a new autosave commit is created.  Subsequent
# changes amend that commit.
#
# MAKE SURE TO RUN THIS FROM A SEPARATE BRANCH FOR EACH INVOCATION!
# Corollary 1:  No two users should be running this on the same branch.
# Corollary 2:  Running on master is not allowed.
#
arg="${1}"

if [ "x${arg}x" = "xx" ] || [ "x${arg}x" = "x-hx" ] || [ "x${arg}x" = "x--helpx" ]; then
    echo "Example usage: "
    echo ""
    echo "    cd local_path_to/top_of_my_git_repo"
    echo "    autosave_continuously.sh ubuntu@devserver.aws.com:remote_path_to/top_of_my_git_repo"
    exit 0
fi

REMOTE_USER=`echo "${arg}"  | awk -F'\(@\|\:\)' '{print $1}'`
REMOTE_HOST=`echo "${arg}"  | awk -F'\(@\|\:\)' '{print $2}'`
REMOTE_WORK_DIR=`echo "${arg}"  | awk -F'\(@\|\:\)' '{print $3}'`

if [ "x`basename ${REMOTE_WORK_DIR}`x" != "x`basename ${PWD}`x" ]; then
    echo "ERROR:  Remote workdir basename `basename ${REMOTE_WORK_DIR}` is different from local workdir basename `basename ${PWD}`."
    echo ""
    $0 --help
    exit 1
fi

echo "Remote target: ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_WORK_DIR}"

branch=`git status | head -1 | sed 's=^On branch =='`
spaceless=`echo $branch | sed 's=[[:blank:]]=='`

if [ "x${branch}x" = "xmasterx" ]; then
    echo "You can't do this on branch 'master'.  Checkout some other branch first."
    exit 1
fi

if [ "x${branch}x" != "x${spaceless}x" ]; then
    echo "Branch name '${branch}' must not contain whitespace."
    exit 2
fi

#git status
#read -p "Press ENTER to start monitoring branch ${branch} for autosaving and autopushing, or CTRL-C to abort." dummy

echo "Autosave/autopush continuous monitoring of local branch ${branch} initiated."
echo "You may press CTRL-C at any time to terminate the monitor."

first_iteration="yes"
autosave_commit_created="no"
sleep=1

while :
do
    /bin/echo -n `date`

    if ! (git status --short | grep '\(^ M \|^A \)') 2>&1 > /dev/null; then
        #
        # No new added or modified files relative to the top commit.
        #
        # First iteration => push to ensure top commit is saved to github,
        # in case edited manually before starting monitor.
        #
        # Subsequent iterations => do nothing unless change detected.
        #
        if [ ${first_iteration} = "no" ]; then
            echo -n -e "\r"
            sleep $sleep
            continue
        else
            echo ""
        fi
    else
        echo ""
        if [ ${autosave_commit_created} = "no" ]; then
            git commit -a --allow-empty --no-edit --quiet -m "Autosave commit created on `date`."
            echo "Created autosave commit on top of branch ${branch}."
            echo "From now on, will amend this commit with your every edit."
            autosave_commit_created="yes"
        else
            git commit -a --amend --allow-empty --no-edit --quiet
            echo "Amended autosave commit on top of branch ${branch} with your changes."
        fi
    fi

    if [ ${first_iteration} = "yes" ]; then
        # This is only useful the first iteration, because that's when the user
        # might have forgotten to add files that won't be saved automatically.
        # TODO:  Track changes to the list of untracked files and warn ONCE
        # when new untracked files appear.
        echo "Note:  Any newly created files will remain untracked and unsaved,"
        echo "       until you CTRL-C and manually git-add them."
        git status --short
    fi

    (git push --force origin ${branch} --quiet ; echo "Pushed to github.") &

    ssh ${REMOTE_USER}@${REMOTE_HOST} "(cd ${REMOTE_WORK_DIR} ; git checkout master --quiet)"

    git push --force --quiet ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_WORK_DIR} ${branch}

    ssh ${REMOTE_USER}@${REMOTE_HOST} "(cd ${REMOTE_WORK_DIR} ;  git checkout ${branch} --quiet)"

    /bin/echo "Pushed to remote host. "

    # If the local directory contains a makefile, do make on the remote host.
    if [ -e Makefile ]; then
        echo "Now building."
        ssh ${REMOTE_USER}@${REMOTE_HOST} "(cd ${REMOTE_WORK_DIR} ; make -j 8 all)" 2>&1 | sed 's=^=REMOTE_HOST:  ='
        echo "Done with remote compilation."
    fi

    wait

    echo "---------------------------------"

    sleep $sleep
    first_iteration="no"

done
