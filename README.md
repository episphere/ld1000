# ld1000

Calculation of linkage distance using the 1000genome project remote data directly. Live at https://episphere.github.io/ld1000 !

## Example usage

### In the local command line

```javascript
await VCF.LDextraction('chr7:16876630','chr7:16863828')
```

### From any location

```javascript
LinkageData = await (await import('https://episphere.github.io/ld1000/export.js')).LDextraction('chr7:16876630','chr7:16863828')
```

Note cache in exportable method VV