var el = document.createElement("div");
el.id = "debian-notice--giac-js-disabled"
el.setAttribute("style", "position: fixed; top: 0; left: 0; \
width: 100%; height: 100%; \
margin: 0; padding: 0; \
background: rgba(128, 128, 128, 0.5); \
z-index: 999; ");
var notice = document.createElement("div");
notice.innerHTML = "The giac.js library is not available in Debian, because emscripten cannot be packaged.<br/>\
See <a href=\"file:///usr/share/doc/giac-doc/README.Debian\">README.Debian</a> in the giac-doc binary package for details.<br/><br/>\
(Click to hide this message.)";
notice.setAttribute("style", "position: relative; \
top: 50%; transform: perspective(1px) translateY(-50%); \
background: white; \
max-width: 34em; \
margin: auto; padding: 0.5em;");
notice.addEventListener("click", function() { el.style.display = "none"; });
el.appendChild(notice);
document.body.insertBefore(el, document.body.childNodes.item(0));
