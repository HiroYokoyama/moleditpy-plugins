/* License metadata shared by the Explorer and its rendering tests. */
function renderLicenseMetadata(value) {
    const label = typeof value === 'string' && value.trim() && value !== 'NOASSERTION'
        ? value.trim() : 'Not specified';
    const escaped = label.replace(/[&<>"']/g, c => ({
        '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
    }[c]));
    return `<div class="meta-row plugin-license">
        <span class="meta-label">
            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
                <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"></path>
                <path d="M14 2v6h6M8 7h2M8 11h8"></path>
                <circle cx="12" cy="16" r="2"></circle>
                <path d="m10.5 17.5-.5 3 2-1 2 1-.5-3"></path>
            </svg>
            License
        </span>
        <span class="meta-value mono" title="Plugin code license. Dependencies may have different terms; see the project notices.">${escaped}</span>
    </div>`;
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { renderLicenseMetadata };
}
