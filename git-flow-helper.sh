#!/bin/bash
# 🚀 gbcms Git Flow Helper Script

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}🚀 gbcms Git Flow Helper${NC}"
echo ""

# Function to create feature branch
create_feature() {
    read -p "Enter feature name (e.g., amazing-new-feature): " feature_name

    if [ -z "$feature_name" ]; then
        echo -e "${RED}❌ Feature name cannot be empty${NC}"
        exit 1
    fi

    branch_name="feature/$feature_name"

    echo -e "${YELLOW}🔄 Creating feature branch: $branch_name${NC}"

    # Ensure we're up to date
    git checkout develop
    git pull origin develop

    # Create and push feature branch
    git checkout -b "$branch_name"
    git push -u origin "$branch_name"

    echo -e "${GREEN}✅ Feature branch created and pushed!${NC}"
    echo -e "${YELLOW}📝 Next: Create a Pull Request to develop branch${NC}"
}

# Function to create release branch
create_release() {
    read -p "Enter version (e.g., 2.1.0): " version

    if [ -z "$version" ]; then
        echo -e "${RED}❌ Version cannot be empty${NC}"
        exit 1
    fi

    branch_name="release/$version"

    echo -e "${YELLOW}🔄 Creating release branch: $branch_name${NC}"

    # Ensure we're up to date
    git checkout develop
    git pull origin develop

    # Create release branch
    git checkout -b "$branch_name"

    echo -e "${GREEN}✅ Release branch created!${NC}"
    echo -e "${YELLOW}📝 Next: Update version numbers and changelog${NC}"
}

# Function to show status
show_status() {
    echo -e "${GREEN}📋 Current Git Status${NC}"
    echo "=========================="
    git status --short
    echo ""
    echo -e "${GREEN}🏷️  Recent commits${NC}"
    git log --oneline -5
}

# Function to cleanup merged branches
cleanup_branches() {
    echo -e "${YELLOW}🧹 Cleaning up merged branches...${NC}"

    # Delete local branches that are merged
    git branch --merged main | grep -v 'main\|develop' | xargs -r git branch -d

    # Show remaining branches
    echo -e "${GREEN}📋 Remaining branches:${NC}"
    git branch -a
}

# Show menu
echo "Available commands:"
echo "1) Create feature branch"
echo "2) Create release branch"
echo "3) Show git status"
echo "4) Cleanup merged branches"
echo "5) Exit"
echo ""

read -p "Choose an option (1-5): " choice

case $choice in
    1)
        create_feature
        ;;
    2)
        create_release
        ;;
    3)
        show_status
        ;;
    4)
        cleanup_branches
        ;;
    5)
        echo -e "${GREEN}👋 Goodbye!${NC}"
        exit 0
        ;;
    *)
        echo -e "${RED}❌ Invalid option${NC}"
        exit 1
        ;;
esac
