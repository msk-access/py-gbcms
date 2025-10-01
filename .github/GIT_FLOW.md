# 🚀 gbcms Git Flow

## 📋 **Branch Structure**

| Branch | Purpose | Protection | Auto-delete |
|--------|---------|------------|-------------|
| **`main`** | Production releases | ✅ Protected | ❌ Never |
| **`develop`** | Development integration | ⚠️ Optional | ❌ Never |
| **`feature/*`** | New features | ❌ None | ✅ On merge |
| **`hotfix/*`** | Production fixes | ❌ None | ✅ On merge |
| **`release/*`** | Release preparation | ❌ None | ✅ On merge |

## 🔄 **Development Workflow**

### **Feature Development**
```bash
# Create feature branch
git checkout develop
git pull origin develop
git checkout -b feature/amazing-new-feature

# Develop your feature
# ... make commits ...

# Push and create PR
git push -u origin feature/amazing-new-feature
# Create Pull Request to develop branch
```

### **Release Process**
```bash
# Create release branch
git checkout develop
git pull origin develop
git checkout -b release/2.1.0

# Update version numbers, changelog
# ... make final adjustments ...

# Merge to main (creates release)
git checkout main
git merge release/2.1.0
git push origin main

# Merge back to develop
git checkout develop
git merge release/2.1.0
git push origin develop

# Delete release branch
git branch -d release/2.1.0
```

### **Hotfix Process**
```bash
# Create hotfix branch from main
git checkout main
git checkout -b hotfix/critical-bug-fix

# Fix the issue
# ... make commits ...

# Merge to main (creates patch release)
git checkout main
git merge hotfix/critical-bug-fix
git push origin main

# Merge to develop
git checkout develop
git merge hotfix/critical-bug-fix
git push origin develop
```

## 🏷️ **Release Workflow Integration**

The **unified release workflow** triggers on:
- **Semantic version tags**: `1.0.0`, `2.1.3`, etc.
- **Manual dispatch** for testing

**Triggered by:**
- Pushing tags matching `[0-9]+.[0-9]+.[0-9]+`
- Manual workflow dispatch

## 📦 **Release Checklist**

### **Before Release**
- [ ] All tests pass
- [ ] Code review completed
- [ ] Documentation updated
- [ ] Version numbers updated
- [ ] Changelog updated

### **Release Steps**
1. **Create release branch** from develop
2. **Final testing** and adjustments
3. **Merge to main** (triggers release workflow)
4. **Tag created** automatically by workflow
5. **PyPI package** published
6. **Docker image** published
7. **Merge back** to develop

## 🔧 **Git Commands Reference**

```bash
# Branch management
git branch -a              # List all branches
git checkout -b feature/x  # Create feature branch
git push -u origin feature/x # Push and track

# Cleanup
git branch -d feature/x    # Delete local branch
git push origin :feature/x # Delete remote branch

# Status
git status                 # Current status
git log --oneline         # Recent commits
git diff                  # Uncommitted changes
```

## 🚨 **Important Notes**

- **Always pull** before creating new branches
- **Use descriptive** branch names with prefixes
- **Create PRs** for all feature branches to develop
- **Test thoroughly** before merging to main
- **Update documentation** with new features

## 🎯 **Next Steps**

1. **Set up branch protection** rules in GitHub
2. **Configure CI** for automatic testing on PRs
3. **Set up release notes** automation
4. **Configure issue templates** for feature requests

**Happy developing!** 🚀
